#
# vim:set fdm=marker:

#################################################################
#Creating Paolas Special tables {{{

#first link the MRM atlas table to the yeast genome {{{
"""
drop table if exists MRMPepLink_final;
create table MRMPepLink_final (
 peptide_key int,
 mrm_key int,
 charge int
);
alter table MRMPepLink_final add index(peptide_key);
"""


c.execute( 'select id,sequence from hroest.MRMAtlas_qtrap_final_no_pyroGlu')
mrm = c.fetchall()
c.execute( 'select id,sequence from ddb.peptide where experiment_key = 3131')
pep = c.fetchall()



#create a mapping of sequence to peptide_id
pepd = {}
for p in pep:
    pepd[p[1]] = p[0]


import re 
res = []
notfound_dict = {}
nontr = 0
other = 0
for m in mrm:
    try:
        spl = m[1].split('/')
        sequence = re.sub('[^A-Z]+', '', m[1] )
        id = pepd[sequence]
        charge = spl[1]
        tuple = ( (m[0], id, charge) )
        res.append( tuple  )
    except: 
        notfound_dict[sequence] = ''
        if sequence.count('K') + sequence.count('R') < 2: print m; other +=1
        else: nontr += 1

len(notfound_dict)
len(res)
#we have around 9511 nontryptic and 2321 other transitions
#we found 323123A => 96.5 %

q = "INSERT INTO hroest.MRMPepLink_final (mrm_key, peptide_key, charge)" + \
     " VALUES (%s,%s,%s) "
c.executemany( q, res)

#27345 peptides
#42599 precursors out of 43822

}}}
#################################################################
#START creating the MRM Atlas tables {{{





exp_key = 3131  #yeast (all)
#exp_key = 3352  #yeast, 1200 peptides
#exp_key = 3445  #mouse (all)

#mrm atlas
all_peptide_query  = """
select distinct peptide.sequence, molecular_weight, ssrcalc,
genome_occurence, 
peptide.id as peptide_key
from peptide 
inner join peptideOrganism on peptide.id = peptideOrganism.peptide_key
inner join compep.ssrcalc_prediction on ssrcalc_prediction.sequence =
peptide.sequence
where experiment_key = %s
and length( peptide.sequence ) > 1
and peptide.sequence in (select sequence from hroest.MRMAtlas_yeast_201009_final)
""" % exp_key





 
 


c2 = db.cursor()
###################################
# A) store the transitions
###################################
#read in the data from DB and put into bins
# 75 peptides/sec
insert_db = True
modify_cysteins = True
read_from_db = False
c.execute( 'use ddb;' )
t = utils.db_table( c2 )
t.read( all_peptide_query )
print "executed the query"
rows = t.c.fetchall()

#1000 entries / 9 s with fast (executemany), 10.5 with index
#1000 entries / 31 s with normal (execute), 34 with index
transition_table = 'hroest.srmTransitions_yeast_mrmatlas'
peptide_table = 'hroest.srmPeptides_yeast_mrmatlas'
c.execute('truncate table ' + peptide_table)
c.execute('truncate table ' + transition_table)
mass_bins = [ []  for i in range(0, 10000) ]
rt_bins = [ []  for i in range(-100, 500) ]
tmp_c  = db.cursor()
tmp2_c  = db.cursor()
progressm = progress.ProgressMeter(total=len(rows), unit='peptides')
start = time.time()
R = silver.Residues.Residues('mono')
exclude = 0
for i,row in enumerate(rows):
    progressm.update(1)
    ttmp = utils.db_table( tmp_c )
    ttmp.read( "select * from hroest.MRMAtlas_yeast_201009_final a where sequence = '%s' and Intensity > 100 order by Intensity " % row[0])
    transitions = ttmp.fetchall_groupBy( 'parent_charge')
    for ch, t_charge in transitions.iteritems(): 
        mod_sequence = ttmp.row( t_charge[0] , 'sequence').split('/')[0]
        #
        peptide = DDB.Peptide()
        peptide.set_sequence( mod_sequence )
        peptide.charge = ttmp.row( t_charge[0] , 'parent_charge')
        peptide.ssr_calc         = t.row( row, 'ssrcalc' )
        peptide.id               = t.row( row, 'peptide_key' )
        #if '[' in peptide.get_modified_sequence('SEQUEST'): continue
        #
        assert ch == peptide.charge
        peptide.create_fragmentation_pattern(R)
        collider.insert_peptide_in_db(peptide, db, peptide_table, transition_group=i)
        prepare  = []
        vals = "type, fragment_number, group_id, q3_charge, q3 "
        q = "insert into %s (%s)" % (transition_table, vals)  + " VALUES (%s,%s,%s,%s,%s)" 
        #attention: we write the parent ion charge into the q3_charge field!
        #this to be able to still figure out the relative order of the
        #transitions later on
        for j,tt in enumerate(t_charge):
            prepare.append( [ttmp.row( tt, 'Ion_description').split('/')[0][:8],
                             j+1, i, ch, ttmp.row( tt, 'q3') ] )
        dummy = tmp2_c.executemany( q, prepare)



#important: also do the isotopes
#A2 create the additional parent ions for the isotope patterns
PATTERNS_UP_TO = 3 #up to how many isotopes should be considered
cursor = c
vals = "peptide_key, q1_charge, q1, modified_sequence, ssrcalc, isotope_nr, transition_group"
query = "SELECT parent_id, %s FROM %s" % (vals, peptide_table)
cursor.execute( query )
allpeptides =  cursor.fetchall()
prepared = []
for p in allpeptides:
    q1_charge = p[2]
    for i in range(1,1+PATTERNS_UP_TO):
        prepared.append( [p[1], p[2], p[3] + (R.mass_diffC13 * i* 1.0) / q1_charge,
                              p[4], p[5], i, p[7]] )

q = "INSERT INTO %s (%s)" % (peptide_table, vals)  + \
     " VALUES (" + "%s," *6 + "%s)"
c.executemany( q, prepared)


#END creating the MRM Atlas tables }}}
#################################################################

#################################################################
#START creating the Peptide Atlas tables {{{
from Bio import SeqIO
f = open('/home/hroest/data/pepatlas/APD_Sc_all.fasta')
peptides = []
for p in SeqIO.parse(f, "fasta"):
    peptides.append(p.seq.tostring())

f.close()

import perl
perl.eval( "use lib '%s'" % perl_ddb_location)
import ddb #needs the perl api to get to the 2DDB
for seq in peptides:
    p = ddb.Peptide(experiment_key = 3452)
    p.peptide = seq 
    p.peptide_type = 'bioinformatics'
    p.addignore_setid()


exp_key = 3452  #yeast (all)
#pep atlas
all_peptide_query  = """
select distinct peptide.sequence, molecular_weight, ssrcalc,
1 as genome_occurence, 
peptide.id as peptide_key
from peptide 
inner join compep.ssrcalc_prediction on ssrcalc_prediction.sequence = peptide.sequence
where experiment_key = %s
and length( peptide.sequence ) > 1
""" % exp_key

c2 = db.cursor()
###################################
# A) store the transitions
###################################
#read in the data from DB and put into bins
# 75 peptides/sec
insert_db = True
modify_cysteins = True
read_from_db = False
c.execute( 'use ddb;' )
t = utils.db_table( c2 )
t.read( all_peptide_query )
print "executed the query"
rows = t.c.fetchall()

#1000 entries / 9 s with fast (executemany), 10.5 with index
#1000 entries / 31 s with normal (execute), 34 with index
transition_table = 'hroest.srmTransitions_yeast_pepatlas_tt'
peptide_table = 'hroest.srmPeptides_yeast_pepatlas_tt'
c.execute('truncate table ' + peptide_table)
c.execute('truncate table ' + transition_table)
mass_bins = [ []  for i in range(0, 10000) ]
rt_bins = [ []  for i in range(-100, 500) ]
tmp_c  = db.cursor()
tmp2_c  = db.cursor()
progressm = progress.ProgressMeter(total=len(rows), unit='peptides')
start = time.time()
#####this si wRONG
R = silver.Residues.Residues('mono')
exclude = 0
transition_group = 0
for row in rows:
    progressm.update(1)
    transition_group += 1 
    #if transition_group >= 1001: break #####FOR TESTING ONLY
    for mycharge in [2,3]:  #precursor charge 2 and 3
        peptide = collider.get_peptide_from_table(t, row)
        peptide.charge = mycharge
        if modify_cysteins: peptide.modify_cysteins()
        peptide.create_fragmentation_pattern(R)
        if insert_db:
            #insert peptide into db
            collider.insert_peptide_in_db(peptide, db, peptide_table,
                                          transition_group=transition_group)
    #we want to insert the fragments only once per peptide
    if insert_db:
        #insert fragment charge 1 and 2 into database
        peptide.mass_H = R.mass_H
        collider.fast_insert_in_db(peptide, db, 1, transition_table,
                                   transition_group=transition_group)
        collider.fast_insert_in_db(peptide, db, 2, transition_table,
                                   transition_group=transition_group)


#A2 create the additional parent ions for the isotope patterns
PATTERNS_UP_TO = 3 #up to how many isotopes should be considered
cursor = c
vals = "peptide_key, q1_charge, q1, modified_sequence, ssrcalc, isotope_nr, transition_group"
query = "SELECT parent_id, %s FROM %s" % (vals, peptide_table)
cursor.execute( query )
allpeptides =  cursor.fetchall()
prepared = []
for p in allpeptides:
    q1_charge = p[2]
    for i in range(1,1+PATTERNS_UP_TO):
        prepared.append( [p[1], p[2], p[3] + (R.mass_diffC13 * i* 1.0) / q1_charge,
                              p[4], p[5], i, p[7]] )



#END creating the Peptide Atlas tables }}}
#################################################################

#################################################################
#START creating the top3 tables {{{
observations, pI, precursor m/z, average retention time, precursor charge, peptide, protein. All peptides are proteotypic.

drop table hroest.hl_top3  ;
create table hroest.hl_top3  (
id int auto_increment primary key, 
observations int, 
total_obs int, 
pI double, 
q1 double, 
rt double, 
q1_charge int, 
sequence varchar(255),
sequence_ch varchar(255),
sequence_naked varchar(255),
protein varchar(255)
)

values = "observations , total_obs , pI , q1 , rt , q1_charge , sequence_ch, sequence , sequence_naked , protein "

import csv
r = csv.reader( open('/home/hroest/data/targets_3pep.csv'))
cursor = c
prepare  = []
for row in r:
    seq = row[5]
    mystrip = "".join( [s for s in seq if s in string.uppercase] )
    obs = row[0].split('/')[0]
    total = row[0].split('/')[1]
    rr = [obs, total]
    rr.extend( row[1:6])
    rr.append(seq.split('/')[0])
    rr.append( mystrip )
    rr.append( row[-1])
    prepare.append( rr )

c.executemany( 
" insert into hroest.hl_top3 (%s) values " % values + 
    ' (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', prepare)


alter table hroest.hl_top3 add index(sequence_naked);
alter table hroest.hl_top3 add index(sequence);
alter table hroest.hl_top3 add index(sequence_ch);

select count(*) from hroest.hl_top3;


#946 are not tryptic
select count(*) from hroest.hl_top3 t inner join ddb.peptide p
on p.sequence = t.sequence_naked where experiment_key = 3131;

#I found 4560
select count(*) from hroest.hl_top3 t inner join hroest.srmPeptides_yeast s
on s.modified_sequence = t.sequence;

#I found 2400
select *
from hroest.hl_top3 t inner join 
hroest.MRMAtlas_yeast_qqq_201002_q2  
a on a.sequence = t.sequence_ch
group by a.sequence


#I found 4560
#create the top3 tables, which only contain all b,y ions 
#for the peptides seen in henry lams top3 estimate
create table hroest.srmPeptides_yeast_top3 as
select 
s.parent_id        ,
s.peptide_key      ,
s.modified_sequence,
s.q1_charge        ,
s.q1               ,
s.ssrcalc          ,
s.isotope_nr       ,
s.transition_group 
from hroest.hl_top3 t inner join hroest.srmPeptides_yeast s
on s.modified_sequence = t.sequence;

alter table hroest.srmPeptides_yeast_top3 add index(peptide_key);
alter table hroest.srmPeptides_yeast_top3 add index(q1);
alter table hroest.srmPeptides_yeast_top3 add index(ssrcalc);
alter table hroest.srmPeptides_yeast_top3 add index(transition_group);

create table hroest.srmTransitions_yeast_top3 as
select
srm_id          ,
group_id        ,
q3_charge       ,
q3              ,
type            ,
fragment_number 
from hroest.srmTransitions_yeast s
where s.group_id  in (select transition_group from hroest.srmPeptides_yeast_top3);
alter table hroest.srmTransitions_yeast_top3 add index(group_id);
alter table hroest.srmTransitions_yeast_top3 add index(q3);

#print "\n".join(chr(i) + " ... : Text" for i in range(65,91,2))

#END creating the top3 tables }}}
#################################################################

#################################################################
#START creating the MRM Atlas (all transitions) tables {{{

#will have matching parent_ids to hroest.srmPeptides_yeast
drop table hroest.srmPeptides_yeast_mrmatlasall ;
create table hroest.srmPeptides_yeast_mrmatlasall as
select  distinct
s.parent_id        ,
s.peptide_key      ,
s.modified_sequence,
s.q1_charge        ,
s.q1               ,
s.ssrcalc          ,
s.isotope_nr       ,
s.transition_group 
from  hroest.MRMPepLink_final t inner join hroest.srmPeptides_yeast s
on s.peptide_key = t.peptide_key;
#from MRMAtlas_yeast_201009_final t inner join hroest.srmPeptides_yeast s
#on t.modified_sequence = CONCAT( s.modified_sequence, '/', s.q1_charge);
alter table hroest.srmPeptides_yeast_mrmatlasall add index(peptide_key);
alter table hroest.srmPeptides_yeast_mrmatlasall add index(q1);
alter table hroest.srmPeptides_yeast_mrmatlasall add index(ssrcalc);
alter table hroest.srmPeptides_yeast_mrmatlasall add index(transition_group);

drop table hroest.srmTransitions_yeast_mrmatlasall ;
create table hroest.srmTransitions_yeast_mrmatlasall as
select
srm_id          ,
group_id        ,
q3_charge       ,
q3              ,
type            ,
fragment_number 
from hroest.srmTransitions_yeast s
where s.group_id  in (select transition_group from hroest.srmPeptides_yeast_mrmatlasall);
alter table hroest.srmTransitions_yeast_mrmatlasall add index(group_id);
alter table hroest.srmTransitions_yeast_mrmatlasall add index(q3);
#END creating the MRM Atlas (all transitions) tables }}}
#################################################################

#################################################################
#{{{ 2010,12,13 creating complete yeast digest + peptide_atlas
import MySQLdb
sys.path.append( '/home/hroest/lib/' )
sys.path.append( '/home/hroest/lib/hlib' )
sys.path.append( '/home/hroest/srm_clashes/code' )
sys.path.append( '/home/hroest/msa/code/tppGhost' ) #DDB
db = MySQLdb.connect(read_default_file="~/.my.cnf")
c = db.cursor()
import silver
import collider
import DDB 
R = silver.Residues.Residues('mono')
import progress

c.execute("select mod_sequence from hroest.peptideatlas_201012_davecampbell")
seqs = c.fetchall()

peptide = DDB.Peptide()
excluded = 0
allowed = ['M[147]', 'A', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M',
               'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y', 'C[160]']
#I excluded 24818 sequences
peptideatlas = []
for s in seqs:
    peptide.set_sequence(s[0])
    seq = s[0]
    excl = False
    fragments =  list(peptide._get_modified_fragments() )
    for fr in fragments:
        if fr not in allowed: excl = True
    if excl: excluded += 1
    else: peptideatlas.append( seq )

peptideatlas = set(peptideatlas)

c.execute('select distinct modified_sequence from hroest.srmPeptides_yeast')
seqs = c.fetchall()
proteome = set( [s[0] for s in seqs])

union = peptideatlas.union(proteome)
len(proteome)
len(peptideatlas)
#I go from 192792 to 222298
#I thus add 29506 peptides
len(union)

transition_table = 'hroest.srmTransitions_yeast_pepatlas_union'
peptide_table = 'hroest.srmPeptides_yeast_pepatlas_union'
c.execute('truncate table ' + peptide_table)
c.execute('truncate table ' + transition_table)
tmp_c  = db.cursor()
progressm = progress.ProgressMeter(total=len(union), unit='peptides')
peptide = DDB.Peptide()
insert_db = True
for ii,sequence in enumerate(union):
    progressm.update(1)
    for mycharge in [2,3]:  #precursor charge 2 and 3
        peptide.set_sequence(sequence)
        peptide.charge = mycharge
        peptide.mass_H = R.mass_H
        peptide.create_fragmentation_pattern(R)
        tmp = c.execute( """
        select p.id from ddb.peptide p 
             where p.sequence = '%s' and experiment_key = 3131
        """ % peptide.sequence
        )
        rr = c.fetchall()
        try: peptide.id = rr[0][0]
        except Exception: peptide.id = -99
        tmp = c.execute( """
        select ssrcalc from compep.ssrcalc_prediction 
             where sequence = '%s' 
        """ % peptide.sequence
        )
        rr = c.fetchall()[0]
        peptide.ssr_calc= rr[0]
        if insert_db:
            #insert peptide into db
            ##collider.insert_peptide_in_db(S, db, peptide_table, transition_group=i)
            transition_group = ii
            c = db.cursor()
            #insert peptide into db
            vals = "peptide_key, q1_charge, q1, ssrcalc, modified_sequence, isotope_nr, transition_group"
            q = "insert into %s (%s) VALUES (%s,%s,%s,%s,'%s', %s, %s)" % (
                peptide_table,
                vals, 
                peptide.id, peptide.charge, 
                peptide.charged_mass,
                peptide.ssr_calc, 
                peptide.get_modified_sequence(),
                0, #we only have the 0th isotope (0 C13 atoms)
                transition_group
            )
            tmp = c.execute(q)
            peptide.parent_id = db.insert_id()
    #we want to insert the fragments only once per peptide
    if insert_db:
        #insert fragment charge 1 and 2 into database
        collider.fast_insert_in_db( peptide, db, 1, transition_table, transition_group=ii)
        collider.fast_insert_in_db( peptide, db, 2, transition_table, transition_group=ii)




#A2 create the additional parent ions for the isotope patterns
PATTERNS_UP_TO = 3 #up to how many isotopes should be considered
cursor = c
vals = "peptide_key, q1_charge, q1, modified_sequence, ssrcalc, isotope_nr, transition_group"
query = "SELECT parent_id, %s FROM %s" % (vals, peptide_table)
cursor.execute( query )
allpeptides =  cursor.fetchall()
prepared = []
for p in allpeptides:
    q1_charge = p[2]
    for i in range(1,1+PATTERNS_UP_TO):
        prepared.append( [p[1], p[2], p[3] + (R.mass_diffC13 * i* 1.0) / q1_charge,
                              p[4], p[5], i, p[7]] )

q = "INSERT INTO %s (%s)" % (peptide_table, vals)  + \
     " VALUES (" + "%s," *6 + "%s)"
c.executemany( q, prepared)



#}}}
#################################################################

create table
hroest.srmPeptides_yeast_pI15
select
parent_id         ,
peptide_key       ,
modified_sequence ,
q1_charge         ,
q1                ,
ssrcalc           ,
isotope_nr        ,
transition_group  
from srmPeptides_yeast p
inner join ddb.peptide d on d.id = p.peptide_key
where pI >= 14 and pI < 15
;


#}}}

drop TABLE hroest.humanpeptides;
CREATE TABLE hroest.humanpeptides (
        id int auto_increment,
        sequence  varchar(255) )
        ;

LOAD DATA LOCAL INFILE
'/tmp/humanpeptides.csv'
INTO TABLE hroest.humanpeptides 
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
#IGNORE 1 LINES 
( sequence )
;


#################################################################

#################################################################
#{{{ MRM Atlas Ion Abundance analysis

c.execute( 
"""
select id, Ion_description 
from hroest.MRMAtlas_yeast_201009_final
"""
)
orig_desc = c.fetchall()

desc = []
for d in orig_desc:
    #we only append the first, e.g. "best" explanation
    dd = d[1].split(',')[0]
    desc.append( [ d[0], dd] )

    #append all explanations => not a good idea
    #for dd in d[1].split(','):
    #    desc.append( [ d[0], dd] )


prepare = []
bra_c = 0
for d in desc:
    id = d[0]
    d = d[1]
    ion = d.split('/')[0]
    #for some reasons, 806 have a [] around them
    # => it might mean that they are wrongly associated
    isotope = 0
    if ion.startswith('['): bra_c += 1; continue
    if ion.find('i') != -1:
        if ion.find('i') != len(ion)-1: break
    if ion.endswith('i'): ion = ion[:-1]; isotope = 1
    charge = ion.split('^')
    if len(charge) == 1: charge = 1
    else: ion = charge[0]; charge = int( charge[1] )
    loss = ion.split('-')
    if len(loss) == 1: loss = 0
    else: ion = loss[0]; loss = int( loss[1] )
    gain = ion.split('+')
    if len(gain) == 1: gain = 0
    else: ion = gain[0]; gain = int( gain[1] )
    fragment = ion[0]
    if len(ion) > 1: nr = int(ion[1:])
    else: nr = 0
    #print  d, fragment, nr, charge, loss, gain, isotope
    prepare.append( [  id, fragment, nr, charge, loss, gain, isotope] )
    #if  id > 20: break


create table hroest.result_mrmfragment_analysis_final_best_expl (
id int, 
fragment char(1), 
nr int,
charge int,
loss int, 
gain int,
isotope int
)

c.executemany( 
   'insert into hroest.result_mrmfragment_analysis_final_best_expl values (%s,%s,%s,%s,%s,%s,%s)', 
              prepare)


desc hroest.MRMAtlas_yeast_qqq_201002_q2

select sequence,Ion_description, 
@ion := LEFT(Ion_description,  INSTR( Ion_description, '/') - 1) as ion,
@fr := LEFT( @ion, 1)  as fr,
@str := LEFT(@ion,  INSTR( @ion, '-') - 1) as str
from hroest.MRMAtlas_yeast_qqq_201002_q2
limit 100;

###there is no fragment a with any gains/losses
select @all := count(*) from hroest.result_mrmfragment_analysis_final_best_expl;
select @precursor := count(distinct modified_sequence) from hroest.MRMAtlas_yeast_201009_final;
drop table tmp ;
create temporary table tmp as
select length(substr(sequence, 1,15)) -3-1  as t
from hroest.MRMAtlas_yeast_201009_final group by modified_sequence ;
select @estimated_nr_fragments := sum(t) from tmp;

select @all := count(*) from hroest.result_mrmfragment_analysis_final_best_expl an
inner join hroest.MRMAtlas_yeast_201009_final at on at.id = an.id
where Intensity > 1000;

select 
#count(*) / @all * 100 as pcnt, 
CONCAT( fragment, '^', charge) as frag,
parent_charge ch, loss, gain , 
count(*)  as occ
from hroest.result_mrmfragment_analysis_final_best_expl an
inner join hroest.MRMAtlas_yeast_201009_final at on at.id = an.id
where Intensity > 2000
group by fragment, parent_charge, charge, loss, gain
order by occ desc
;

select 
count(*)
from hroest.result_mrmfragment_analysis_final an
#inner join hroest.MRMAtlas_yeast_201009_final at on at.id = an.id
;



c = cursor
c.execute ( """
select 
#count(*) / @all * 100 as pcnt, 
CONCAT( fragment, '^', charge) as frag,
parent_charge ch, loss, gain , 
count(*)  as occ
from hroest.result_mrmfragment_analysis_final_best_expl an
inner join hroest.MRMAtlas_yeast_201009_final at on at.id = an.id
where Intensity > 2000
group by fragment, parent_charge, charge, loss, gain
order by occ desc
; """ )
result = c.fetchall()

import latex
total = sum( [r[4] for r in result])
newr = []
cumm = 0
for r in result:
    rr = list(r)
    cumm += r[4]
    rr.extend( [cumm, cumm*1.0/total])
    newr.append(rr)

tbl = ''
for r in newr:
    tbl += ' ' + latex.get_latex_row([ r'\verb|' + r[0] + '|' , r[1], r[2], r[3], 
                              latex.prepare_int( r[4], 6) ,
                              latex.prepare_int( r[5], 6) ])   + '\n'

mytable = latex.get_longtable(tbl, len(r), "Ion & Precursor Charge & Neutral Loss & Neutral Gain & Occurence & Cummulative Occurence", 
          "Fragment Ion Abundances.", 
          """ Abundance of individual fragment ions as determined on a synthetic 
          set of 28357 peptides on a QQQ machine.
          """, 'table:ion_abundance')

f = open( '/home/hroest/phd/papers/srm_clashes/table_ion_abundance_data.tex', 'w')
f.write( tbl ) 
f.close()
f = open( '/home/hroest/phd/papers/srm_clashes/table_ion_abundance.tex', 'w')
f.write( mytable ) 
f.close()


}}}
{{{#some analysis of the peptide atlas data

#peptide Atlas
#has key 3452
cursor.execute(
"""
create table hroest.delete as select *  from ddb.peptide where experiment_key = 3452;
create table hroest.patlasNoGenome as select * from hroest.delete d where
    d.sequence not in (select sequence from ddb.peptide where experiment_key =
    3131);
"""
)
cursor.execute("select sequence from hroest.patlasNoGenome")
list = cursor.fetchall()
cursor.execute("select sequence from hroest.patlasNoGenome")

nontr = 0
other = []
for s in list:
    seq = s[0]
    if seq.count('K') + seq.count('R') > 1: nontr += 1
    else: other.append( seq)
    

#total 25075
select count(*) from hroest.patlasNoGenome; 
#22139 are tryptic
select count(*) from hroest.patlasNoGenome where
sequence like '%K' or
sequence like '%R' ;
#6540 contain no missed cleavages at the end
drop table hroest.patlasNoGenome2 ;
create table hroest.patlasNoGenome2 as select * 
#select count(*) 
from hroest.patlasNoGenome where
(sequence like '%K' or
sequence like '%R' ) and not  (
sequence like '%K%K' or
sequence like '%K%R' or
sequence like '%R%R' or
sequence like '%R%K');

cursor.execute("select sequence from hroest.patlasNoGenome2;")
nf = cursor.fetchall()
cursor.execute("select sequence from ddb.peptide where experiment_key = 3131")
all = cursor.fetchall()

count = 0
others = []
for n in nf:
    ishere = False
    for a in all:
        if a[0].find(n[0]) != -1: print a, n; ishere = True
    if ishere: count += 1
    else: others.append( n[0] )


len(nf)
count

DIENQYETQITQIEHEVSSSGQEVQSSAK
cursor.execute("drop table if exists hroest.patlasNoGenome3" )
cursor.execute("create table hroest.patlasNoGenome3 (sequence varchar(255) )" )
cursor.executemany( 'insert into hroest.patlasNoGenome3 (sequence) values (%s)', others)

#225 are human contaminants
#select count(*) from ddb.peptide 
drop table hroest.patlasNoGenome4 ;
create table hroest.patlasNoGenome4 as
select sequence from hroest.patlasNoGenome3 pa
where pa.sequence not in (select sequence from ddb.peptide where experiment_key = 3130)




select * from hroest.patlasNoGenome limit 300;






}}}
#################################################################

#################################################################
{{{ #Analysis of the Manndata


manndata_lfq : normalized version
manndata_lfq_sorted : nonnormalized version
tppsearchresults : TPP results of manndata for spectral counting

#make map table
drop table if exists tmp_pepmapping;
create temporary table tmp_pepmapping as
select distinct peptide_key, modified_sequence from 
hroest.srmPeptides_yeast ;
alter table tmp_pepmapping add index(peptide_key);
alter table tmp_pepmapping add index(modified_sequence);


#desc manndata_lfq;
#select count(*) from manndata_lfq;
drop table if exists tmp_peptide ; 
create temporary table tmp_peptide as 
select peptide, sum(obs) as nr, max(average) as av , count(*) as nr_triplicates
from manndata_lfq_sorted
   where obs != 0
   group by peptide;



#in total 713'715 features mapped to identifications
select sum(obs) from manndata_lfq;
#we have 33k peptides
select count(*) from tmp_peptide;
#we have 26k tryptic peptides
select count(*) from tmp_peptide t 
inner join hroest.tmp_pepmapping p on
t.peptide = p.modified_sequence ;
#24k proteotypic peptides
select count(*) from tmp_peptide t 
inner join hroest.tmp_pepmapping p on t.peptide = p.modified_sequence 
inner join ddb.peptideOrganism o on o.peptide_key = p.peptide_key
where genome_occurence = 1
;


#we have 3687 unique proteins
#with maxquant they have 

drop table tmp_allfoundproteins ;
create temporary table tmp_allfoundproteins as 
select count(*) as nr_peptides, protein_key from  tmp_peptide t 
inner join tmp_pepmapping map on map.modified_sequence = t.peptide 
inner join ddb.peptide p on map.peptide_key = p.id
inner join ddb.peptideOrganism o on o.peptide_key = p.id
inner join ddb.protPepLink l on l.peptide_key = p.id
where experiment_key = 3131
and genome_occurence = 1
group by protein_key
order by count(*);

#we have 2305 proteins of which we have 3 or more tryptic peptides
select count(*) from tmp_allfoundproteins where nr_peptides > 2;
select count(*) from tmp_allfoundproteins where nr_peptides = 2;
select count(*) from tmp_allfoundproteins where nr_peptides = 1;

#estimate FDR
select count(*) *0.01*0.01*0.01 from tmp_allfoundproteins where nr_peptides > 2;
select count(*) *0.01*0.01 from tmp_allfoundproteins where nr_peptides = 2;
select count(*) *0.01 from tmp_allfoundproteins where nr_peptides = 1;

#from 27k peptides, we have 8k in the lab
drop table if exists tmp_mannanalysis_pepperprot ;
create table tmp_mannanalysis_pepperprot as
select t.peptide, av as intensity, genome_occurence, protein_key from  tmp_peptide t 
inner join tmp_pepmapping map on map.modified_sequence = t.peptide 
inner join ddb.peptide p on map.peptide_key = p.id
inner join ddb.peptideOrganism o on o.peptide_key = p.id
inner join ddb.protPepLink l on l.peptide_key = p.id
#inner join lab_peptides lab on p.sequence = lab.sequence
where experiment_key = 3131
and genome_occurence = 1
order by protein_key, intensity
;




cursor.execute( 'select * from hroest.tmp_mannanalysis_pepperprot')
mypep = cursor.fetchall()

cursor.execute( 'select * from hroest.lab_peptides')
lab = cursor.fetchall()

labs = [l[3] for l in lab]

proteins = {}
for p in mypep:
    if proteins.has_key(p[3]): proteins[p[3]].append( p)
    else: proteins[p[3]] = [p]


for nr in range(1,7):
    intop = []
    for k,v in proteins.iteritems():
        intop.extend( v[-nr:])
    #we have 9234 of the top3 peptides
    #679 + 2*469 + 2539*3
    haveinlab = 0
    for p in intop:
        if p[0] in labs: haveinlab += 1
    print nr, haveinlab, len(intop), haveinlab*1.0 / len(intop)

#of those 9234 we have 4393

#of the top1 (3687) we have 1903 = 51%
#of the top2 (6695) we have 3326 = 50%
#of the top3 (9234) we have 4393 = 47%
#of the top4 (11420) we have 5186 = 45%
#of the top5 (13313) we have 5821 = 43%
#of the top6 (14965) we have 5276 = 41%


#TODO run mayo on my data:
    go to peptides.csv, get all proteotypic peptides (only one n_proteins)
    then write a file with 
    specid, peptide, protid, modification, score
    dummy,  PEPTIDE, DECOY_da, none, dummy

=> then estimate FDR of the protein identifications


#{{{tinas proteins
#select sequence_key, ac from ddb.isbAc where ac in ( 
select ac, t.peptide, 
t.av as max_triplicate_log_intensity,
t.nr as nr_files_observed, 
genome_occurence
from  tmp_peptide t 
inner join tmp_pepmapping map on map.modified_sequence = t.peptide 
inner join ddb.peptide p on map.peptide_key = p.id
inner join ddb.peptideOrganism o on o.peptide_key = p.id
inner join ddb.protPepLink l on l.peptide_key = p.id
inner join ddb.protein prot on l.protein_key = prot.id
inner join ddb.isbAc ac on ac.sequence_key = prot.sequence_key
where p.experiment_key = 3131
and ac.ac in 
(
'YBR249C',
'YCL017C',
'YDR283C',
'YEL031W',
'YER178W',
'YFL030W',
'YGL006W',
'YGL248W',
'YGL253W',
'YGR193C',
'YGR240C',
'YHR107C',
'YHR166C',
'YHR183W',
'YIL084C',
'YJL026W',
'YJL136C',
'YJR051W',
'YKL060C',
'YKL145W',
'YKR031C',
'YLR039C',
'YLR058C',
'YLR249W',
'YLR403W',
'YML109W',
'YMR037C',
'YMR170C',
'YMR205C',
'YMR220W',
'YNR067C',
'YOL116W',
'YPR118W')
group by ac,peptide
order by ac, av DESC
;



#now check the peptide atlas
select p.sequence_key, ac, pep.sequence from ddb.isbAc  ac
inner join ddb.protein p  on p.sequence_key = ac.sequence_key
inner join ddb.protPepLink l on l.protein_key = p.id
inner join ddb.peptide pep on l.peptide_key = pep.id
inner join ddb.peptideOrganism o on o.peptide_key = pep.id
where p.experiment_key = 3131
and genome_occurence = 1
and  pep.sequence in (select sequence from ddb.peptide where experiment_key = 3452)
and ac in ( 
'YBR249C',
'YCL017C',
'YDR283C',
'YEL031W',
'YER178W',
'YFL030W',
'YGL006W',
'YGL248W',
'YGL253W',
'YGR193C',
'YGR240C',
'YHR107C',
'YHR166C',
'YHR183W',
'YIL084C',
'YJL026W',
'YJL136C',
'YJR051W',
'YKL060C',
'YKL145W',
'YKR031C',
'YLR039C',
'YLR058C',
'YLR249W',
'YLR403W',
'YML109W',
'YMR037C',
'YMR170C',
'YMR205C',
'YMR220W',
'YNR067C',
'YOL116W',
'YPR118W')
group by ac, pep.sequence
order by ac
;

















ac, t.peptide, 
t.av as max_triplicate_log_intensity,
t.nr as nr_files_observed, 
genome_occurence


{{{ new analysis Feb 16th

drop table tmp_tinasproteins ;
create table tmp_tinasproteins  as
#tinas proteins
#select sequence_key, ac from ddb.isbAc where ac in ( 
select  p.id as peptide_key, p.sequence as sequence, ac.ac as protein , 
o.genome_occurence as go
from  ddb.peptide p 
inner join ddb.peptideOrganism o on o.peptide_key = p.id
inner join ddb.protPepLink l on l.peptide_key = p.id
inner join ddb.protein prot on l.protein_key = prot.id
inner join ddb.isbAc ac on ac.sequence_key = prot.sequence_key
where p.experiment_key = 3131
and genome_occurence = 1
and ac.ac in 
( 
'YBR249C',
'YCL017C',
'YDR283C',
'YEL031W',
'YER178W',
'YFL030W',
'YGL006W',
'YGL248W',
'YGL253W',
'YGR193C',
'YGR240C',
'YHR107C',
'YHR166C',
'YHR183W',
'YIL084C',
'YJL026W',
'YJL136C',
'YJR051W',
'YKL060C',
'YKL145W',
'YKR031C',
'YLR039C',
'YLR058C',
'YLR249W',
'YLR403W',
'YML109W',
'YMR037C',
'YMR170C',
'YMR205C',
'YMR220W',
'YNR067C',
'YOL116W',
'YPR118W'
)
order by ac DESC
;


desc tmp_tinasproteins ;



#tinas proteins
#select sequence_key, ac from ddb.isbAc where ac in ( 
select tina.protein, t.peptide, tina.peptide_key,
t.av as max_triplicate_log_intensity,
t.nr as nr_files, t.nr_triplicates ,
tina.go
from  tmp_peptide t 
inner join tmp_pepmapping map on map.modified_sequence = t.peptide 
inner join tmp_tinasproteins tina on tina.peptide_key = map.peptide_key
order by protein
;

desc manndata_lfq;
select * from manndata_lfq limit 1;
desc tmp_pepmapping;

select * from manndata_lfq 
where peptide = 'HNPTSDVATLYSLVNWESTK'
and obs != 0
;


select * from manndata_lfq 
where peptide in
(select sequence from tmp_tinasproteins where protein ='YNR067C')
and obs != 0
;



select sequence from tmp_tinasproteins where protein = 'YML109W';







cursor.execute( """use hroest""")
cursor.execute( """ select distinct protein from tmp_tinasproteins ; """)
import csv
file = open('/tmp/cludwig.csv', 'w')
w = csv.writer(file)
pfile = open('/tmp/cludwig_order.csv' , 'w')
pw = csv.writer(pfile)
pw.writerow( ['Protein', 'Sequence', 'Average', 'Std Deviation', 'Nr observations', 'File'])
r = ['Protein', 'MaxOcc', '']
r.extend( range(1,76) )
w.writerow(r)
r = ['','','']
for i in range(5): r.extend( range(1,16))

w.writerow(r)
for protein in cursor.fetchall():
    cursor.execute(
    """
    select sequence from tmp_tinasproteins where protein ='%s'
    """ % protein[0] )
    prepareseq = ''
    for peptide in cursor.fetchall():
        prepareseq += "'%s',\n" % peptide[0]
    prepareseq = prepareseq[:-2]
    #
    q = """
    select count(*), file
    from manndata_lfq_sorted
    where peptide in
    (%s)
    and obs != 0
    group by file
    ;
    """ % prepareseq
    cursor.execute( q )
    tmpres = cursor.fetchall()
    r = [protein[0]]
    mymax = 0
    if len(tmpres) > 0:mymax= max( [res[0] for res in tmpres])
    r.append(mymax)
    r.append('')
    mymap = [0 for i in range(75)] 
    for res in tmpres:
        #r.append( res[0] )
        mymap[ res[1] ] = res[0]
    r.extend(mymap)
    w.writerow(r)
    #
    q = """
    select '%s', peptide, average, stdev, obs, file
    from manndata_lfq_sorted
    where peptide in
    (%s)
    and obs != 0
    order by file,average DESC
    ;
    """ % (protein[0], prepareseq)
    cursor.execute(q)
    pw.writerows( cursor.fetchall() )

pfile.close()
file.close()



cursor.execute( """use hroest""")
cursor.execute( """ select distinct protein from tmp_tinasproteins ; """)
import csv
file = open('/tmp/cludwig_speccount.csv', 'w')
w = csv.writer(file)
r = ['Protein', 'MaxOcc', '']
r.extend( range(1,76) )
w.writerow(r)
w.writerow([])
r = ['','','']
for i in range(5): r.extend( range(1,16))

w.writerow(r)
for protein in cursor.fetchall():
    cursor.execute(
    """
    select sequence from tmp_tinasproteins where protein ='%s'
    """ % protein[0] )
    prepareseq = ''
    for peptide in cursor.fetchall():
        prepareseq += "'%s',\n" % peptide[0]
    prepareseq = prepareseq[:-2]
    #
    q = """
    select count(*), file
    from tmp_tppsearch_analysis_2
    where sequence in
    (%s)
    group by fraction, sds
    ;
    """ % prepareseq
    cursor.execute( q )
    tmpres = cursor.fetchall()
    r = [protein[0]]
    mymax = 0
    if len(tmpres) > 0:mymax= max( [res[0] for res in tmpres])
    r.append(mymax)
    r.append('')
    mymap = [0 for i in range(75)] 
    for res in tmpres:
        #r.append( res[0] )
        mymap[ res[1] ] = res[0]
    r.extend(mymap)
    w.writerow(r)





    q = """
    select peptide, average, stdev, obs
    from manndata_lfq 
    where peptide in
    (%s)
    and obs != 0
    and file = 28
    order by file, average DESC
    ;
    """ % prepareseq
    cursor.execute(q)
    newr = cursor.fetchall()


this = newr[0]
next = newr[1]

for i in range(1,len(newr)):

for this,next in zip(newr[:-1], newr[1:]):
    thisvar = this[2] * this[2]
    nextvar = next[2] * next[2]
    combined_s2 =  thisvar*this[3] + nextvar*next[3] / (this[3] + next[3] - 2)
    ttest = (this[1] - next[1]) / numpy.sqrt(combined_s2 * (1.0/this[3] + 1.0/next[3]) )
    print this[0], next[0], ttest

for this,next in zip(newr[:-1], newr[2:]):
    thisvar = this[2] * this[2]
    nextvar = next[2] * next[2]
    combined_s2 =  thisvar*this[3] + nextvar*next[3] / (this[3] + next[3] - 2)
    ttest = (this[1] - next[1]) / numpy.sqrt(combined_s2 * (1.0/this[3] + 1.0/next[3]) )
    print this[0], next[0], ttest



select count(*), file
from manndata_lfq 
where peptide in
(
'MQLYLTLLFLLSFVECSYISFISNNADEILETDLIETLSYATLTVGEPYVAQSVVVTR',
'VSAASHSPLSVSPK',
'VSASPINSQDSDSNTR',
'TAVQLSLSLSNYASQVSQK',
'ISAQTNNDPVTVSNIYANDNSK',
'SSVHNLSSVSGVASVMPSASTMR',
'VTTLLSQTASTSTSTLFSSSLSISGTQLNGTLLTSVSK',
'GTIDPLVTQMPSYSSQETK',
'IIPSSLTSNK',
'TIYTISVR',
'TNAATATGEDSFIASTPASSTLFYPSNSTQDLVQTLASTTASPAYPSNR',
'TQITLSPSVSLYSTTSPIYPSNITENGSSPSPSLSSTVSPVYPSSSTGNILLSSLFSTVDSSSSPVSSTLDTIYVSSSMQATISSSSSSR',
'TSSSSLSTSTSSTATTTENSSTTTIVNLFNAVSTDEPPTVFDR',
'SPNPMSLADGVSNDGPIQTNK',
'FYTNLIVGSQESPAFVYPYSLWK',
'YTSSSYGFAVQHTTVDQYSYGGYDSSGNAEYLVNPLGIAHVVFSASNFDSSMTMQVDEMTLSSTR',
'VVLSESNDSSNYLEIPLVQGMGFATGIYHGSLNAK',
'IGSSVGFNTIVSESSSNLAQGILK',
'ITLLNGVTWLCYVIGPDDLTSTDFSLEVSSEYEIK',
'ASASVDGLIIQLAVAPSETDYEVFYDQAAGMYVTNFK',
'LQGVSDGSTATYEFSYTTQGESASGSTMIFALPHHESSFSDIMQDYYTGIQLASTTK',
'GVMNGYLTTSLQFSTSLNR',
'QISWLPWSSQLGSNLLEYSK',
'EQLQLLAEVANSELQVSISESISGLNTYYLGK',
'YSYILLTVSEIIQDEASTK',
'STLENIK',
'SAFDILLQNEQTYPLIYDTK',
'FNGLVSSGDWGSTSTQYDFGNTYYNDHHFHYGYIIHAAAVIGYVDSK',
'LNGTWAADNK',
'DWVNSLVR',
'DVANPSEK',
'DEYFAQSR',
'MFDWFNGHSWAAGLYENGNGK',
'NEESSSEDYNFAYAMK',
'LWGATIGDQSMELR',
'GDLMISIMK',
'DAMNDYFYYQNDNTVEPEEIIGNK',
'VSGILFDNIIDYTTYFGTNTEYIHGIHMLPITPVSSNIR',
'SETFVEEEWQTK',
'IEPIIESIESGWTGILK',
'LNQALFDPVDSYAFFSDSTFDSSTYLDNGMSR',
'TWALAFSGGLANSIA'
)
and obs != 0
group by file
;



#now go and verify

YPR118W


}}}



{{{Ruedis idea for a plot: nr peptides / mol weight in Da for each fraction
  

desc tmp_tppsearch_analysis ;

   #this is spectral counting, otherwise group by sequence, sds, fraction
select sds, fraction, t.sequence, prot.id as pr_key,
   m.sequence_key, m.avg_mass
from hroest.tmp_tppsearch_analysis t 
inner join ddb.peptide p on t.sequence = p.sequence
inner join ddb.protPepLink l on l.peptide_key = p.id
inner join ddb.protein prot on l.protein_key = prot.id
inner join hroest.prot_molmass m on m.sequence_key = prot.sequence_key
#inner join ddb.isbAc ac on ac.sequence_key = prot.sequence_key
#
where sds = 4
and fraction =  'N80'
and p.experiment_key = 3131
limit 100;
  

   {{{ prot_molmass

   create table hroest.prot_molmass (
   sequence_key int primary key,
   avg_mass double
   )

   #first we need the mol. weights of the proteins
   cursor.execute(
       """
       select pr.id, sequence, sequence.id from ddb.protein pr
       inner join ddbMeta.sequence on sequence.id = pr.sequence_key
       where experiment_key = 3131
       """
   )

   alter table hroest.prot_molmass add index(sequence_key)

   proteins = cursor.fetchall()

   seq = proteins[0]

   for seq in proteins:
       #R = silver.Residues.Residues('mono')
       R = silver.Residues.Residues('average')
       import DDB, silver
       p = DDB.Peptide()
       p.set_sequence( seq[1])
       p.charge=1
       p.create_fragmentation_pattern(R)
       cursor.execute(
           """
           insert into hroest.prot_molmass (sequence_key, avg_mass) VALUES
           (%s,%s)
           """ % (seq[2], p.molecular_weight)
       )
   }}}


resdir  = '/home/hroest/srm_clashes/results/manndata/r200lb'
stats = dict([ [fr,[]] for fr in ['MEM', 'N80', 'N60', 'NSol', 'PE' ]])
for fr in ['MEM', 'N80', 'N60', 'NSol', 'PE' ]:
    for sds in range(1,16):
        cursor.execute(
           """
           #this is spectral counting, otherwise group by sequence, sds, fraction
        select 
           #sds, fraction, t.sequence, prot.id as pr_key, m.sequence_key, 
            m.avg_mass / 1000 as mass_da
        from hroest.tmp_tppsearch_analysis t 
        inner join ddb.peptide p on t.sequence = p.sequence
        inner join ddb.peptideOrganism o on o.peptide_key = p.id 
        inner join ddb.protPepLink l on l.peptide_key = p.id
        inner join ddb.protein prot on l.protein_key = prot.id
        inner join hroest.prot_molmass m on m.sequence_key = prot.sequence_key
        #inner join ddb.isbAc ac on ac.sequence_key = prot.sequence_key
        #
        where sds = %s
        and fraction =  '%s'
        and o.genome_occurence = 1
        and p.experiment_key = 3131
           """ % (sds, fr)
        )
        thisexp = cursor.fetchall()
        mmass = [t[0] for t in thisexp]
        stats[fr].append( ( numpy.average(mmass), numpy.std(mmass) ) )
        h, n = numpy.histogram( mmass, 100, (0,200))
        h[-1]  = 0
        gnuplot.Gnuplot.draw_boxes_from_data( [h,n], 
                    output= resdir + 'mwplot_%s_%s.eps' % (fr, sds) ,
                    xlabel = '', ylabel = '', tmp_csv='/tmp/gcsv.csv', keep_data = True)



"""
Correspodning GNUPLOT script
1 set terminal postscript color enhanced solid
2 set border 31 linewidth .3
3 set xlabel ""
4 set ylabel ""
5 set output "/home/hroest/srm_clashes/results/manndata/r200lbmwplot_PE_15.eps"
6 plot "/tmp/gcsv.csv" notitle with boxes lt -1 lw 2
"""


width = 0.21
latex = ''
for fr in ['MEM', 'N80', 'N60', 'NSol', 'PE' ]:
    latex += " \\begin{figure} \n"
    for sds in range(1,16):
        latex += "\subfloat[][Fraction %s (%0.2f $\\pm$ %0.2f)]{\includegraphics[width=%s\\textwidth, angle=-90]{%s}}\n" % (
            sds, 
            stats[fr][sds-1][0], stats[fr][sds-1][1] ,
            width, 
            resdir + 'mwplot_%s_%s.pdf' % (fr, sds), 
            )
        if sds % 3 == 0: latex = latex[:-1] + '\\\\\n'
    latex += "\caption{Figure for %s.}\n"  % fr
    latex += "\end{figure}\n\n"


f = open('/tmp/latex', 'w')
f.write(latex)
f.close()


#which proteins are uniquely found in fraction 1 of the PE fraction?


#this is spectral counting, otherwise group by sequence, sds, fraction
drop table tmpfor1 ; create table tmpfor1 as
select sds, fraction, t.sequence, prot.id as pr_key, m.sequence_key,
m.avg_mass / 1000 as mass_da 
from hroest.tmp_tppsearch_analysis t 
inner join ddb.peptide p on t.sequence = p.sequence 
inner join ddb.peptideOrganism o on o.peptide_key = p.id 
inner join ddb.protPepLink l on l.peptide_key = p.id 
inner join ddb.protein prot on l.protein_key = prot.id 
inner join hroest.prot_molmass m on m.sequence_key = prot.sequence_key
where sds = 1 and fraction =  'PE' and p.experiment_key = 3131 
#and m.avg_mass < 80000
and o.genome_occurence = 1
#and prot.id not in (select pr_key from tmp_pe_notlow3) 
;
select count(*), avg(mass_da), std(mass_da) from tmpfor1;




alter table tmp_pe_notlow3 add index(pr_key);

#this is spectral counting, otherwise group by sequence, sds, fraction
create table tmp_pe_notlow3 as
select prot.id as pr_key
from hroest.tmp_tppsearch_analysis t 
inner join ddb.peptide p on t.sequence = p.sequence
inner join ddb.protPepLink l on l.peptide_key = p.id
inner join ddb.protein prot on l.protein_key = prot.id
inner join hroest.prot_molmass m on m.sequence_key = prot.sequence_key
#inner join ddb.isbAc ac on ac.sequence_key = prot.sequence_key
#
where sds < 4
and fraction =  'PE'
and p.experiment_key = 3131
;



}}}


{{{ 2. Spectral count


#file spectralcount.txt
select count(*) spectralcount,  tina.protein, tina.sequence as sequence, 
tina.peptide_key , genome_occurence as genome_occurnce, 
case when l.mrm_key is NULL then 'no' else 'yes' end
#ifnull(l.mrm_key, 'no'), 
#l.mrm_key as present_in_SRMAltas
#q3, Intensity, Ion_description, parent_charge, m.sequence as modseq, collision_energy,
#m.protein_name as mrmatlas_protein, tina.protein as protein
from hroest.tmp_tinasproteins tina
inner join hroest.tppsearchresults mann on mann.sequence = tina.sequence
inner join ddb.peptideOrganism po on po.peptide_key = tina.peptide_key
left join hroest.MRMPepLink_final l on l.peptide_key = tina.peptide_key
group by tina.sequence
having genome_occurence = 1
order by tina.protein, spectralcount DESC;






    }}}
{{{ 3. Best triplicate and 4. topological sort

cursor.execute('select id, file from hroest.manndata_lfq_filemapping ')
lfqmap = dict(cursor.fetchall() )
cursor.execute('select file, id from hroest.manndata_lfq_filemapping ')
lfqmaprev = dict(cursor.fetchall() )



class Minimal: pass

self = Minimal()
import csv
file = open('/tmp/cludwig_besttriplicate.csv', 'w')
w = csv.writer(file)
w.writerow( ['Protein', 'Sequence', 'Average', 'Std Deviation',   'Nr observations', 'File'])
file2 = open('/tmp/cludwig_topsort.csv', 'w')
w2 = csv.writer(file2)
w2.writerow( ['Protein', 'Sequence',  'Student_tvalue_tonext'])
cursor.execute( """ select distinct protein from hroest.tmp_tinasproteins ; """)
for protein in cursor.fetchall():
    tmp = cursor.execute(
    """
    select sequence from hroest.tmp_tinasproteins where protein ='%s'
    """ % protein[0] )
    prepareseq = ''
    for peptide in cursor.fetchall():
        prepareseq += "'%s',\n" % peptide[0]
    prepareseq = prepareseq[:-2]
    #
    q = """
    select '%s', sequence,name,  sds, fraction, count(*) as obs,
    #left(common_name, length(common_name)-2) as replicate
    LEFT(common_name, length(common_name) - locate("-", reverse(common_name)) ) replicate
    from hroest.tmp_tppsearch_analysis 
    where sequence in
    (%s)
    group by sequence, fraction, sds;
    """ % (protein[0], prepareseq)
    tmp = cursor.execute(q)
    res = cursor.fetchall()
    tppcount = [ [] for i in range(75) ]
    for r in res:
        tppcount[ lfqmaprev[ r[-1] ] ].append(r)
    #
    #
    #
    q = """
    select '%s', peptide, average, stdev,
    #case when l.mrm_key is NULL then 'no' else 'yes' end,
    obs, file
    from hroest.manndata_lfq_sorted mann
    #inner join hroest.tmp_tinasproteins tina on tina.sequence = mann.peptide
    #left join hroest.MRMPepLink_final l on l.peptide_key = tina.peptide_key
    where peptide in
    (%s)
    and obs != 0
    order by file,average DESC
    ;
    """ % (protein[0], prepareseq)
    tmp = cursor.execute(q)
    res = cursor.fetchall()
    #only use those peptides to quantify that were also identified in the TPP
    #search in the same triplicate
    files  = [ [] for i in range(75) ]
    for r in res:
        if r[1] in [s[1] for s in tppcount[ r[-1] ]]:
            files[ r[-1] ].append(r)
    #[len(f) for f in files]
    #[len(f) for f in tppcount]
    #
    mymax = max([len(f) for f in files])
    maxrepl = [f for f in files if len(f) == mymax][0]
    print "number of max", len( ([i for i,f in enumerate(files) if len(f) == mymax]) )
    if len( ([i for i,f in enumerate(files) if len(f) == mymax]) ) >1: 
        print [len(f) for i,f in enumerate(files) ] 
    maxrepl.sort( lambda x,y: -cmp(x[2], y[2] ) )
    w.writerows(maxrepl)
    #
    #
    self.files = files
    sorted = get_sorted(self)
    for a,b in zip(reversed(sorted[:-1]), reversed(sorted[1:])):
        evidence = [e[2] for e in self.edges if e[0] == a and e[1] == b ]
        if len(evidence) > 0: w2.writerow( [protein[0], self.noderevdict[b], max(evidence)] )
        else: w2.writerow( [protein[0], self.noderevdict[b], 0] )
    if len(sorted) > 0: w2.writerow( [protein[0], self.noderevdict[a], 0] )


file.close()
file2.close()

def get_sorted(self):
    files = self.files

    nodes = []
    for f in files:
        for node in f:
            if not node[1] in nodes: nodes.append( node[1] )
    nodedict = dict([(f,i) for i,f in enumerate(nodes) ] )
    noderevdict = dict([(i,f) for i,f in enumerate(nodes) ] )
    edges = []
    for f in files:
        for i,first in enumerate(f):
            for j,second in enumerate(f[:i]):
                #calculate degrees of freedom
                var1 = first[3]**2
                var2 = second[3]**2
                obs1 = first[4]
                obs2 = second[4]
                if obs1 < 2 or obs2 <2: continue
                df = (var1/obs1 + var2/obs2 ) **2
                df /= (var1/obs1)**2/(obs1-1.0) + (var2/obs2)**2/(obs2-1.0)
                tval = (first[2] - second[2] ) / numpy.sqrt( var1/obs1 + var2/obs2 )
                #
                # 100-scipy.stats.t.cdf( tval, df ) * 100
                #
                edges.append(  (nodedict[first[1]], nodedict[second[1]], -tval ) )
                #if tval < -100: print tval, j,i, first, second
    import topsort
    edges.sort( lambda x,y : -cmp(x[2], y[2] ) )
    excluded = []
    runagain = True
    sorted  = []
    while runagain:
        runagain = False
        try:
            myedges = [ (e[0], e[1]) for e in edges if
                        (e[0], e[1]) not in excluded ]
            for i in range(len(myedges)):
                sorted = topsort.topsort( myedges[:i])
        except topsort.CycleError:
            print "excluded" ,i, "/", len(edges),  [e for e in edges if (e[0], e[1]) == myedges[i-1] ]
            excluded.append( myedges[i-1] )
            if i < 1000: runagain = True
    self.edges = edges
    self.nodedict = nodedict
    self.noderevdict = noderevdict
    return sorted





drop table ttmp;
create temporary table  ttmp as
select  name, tina.protein, tina.sequence  , sds, fraction #, fmap.id file
, @replicate := left(common_name, length(common_name)-2) as replicate, 
RIGHT(@replicate, 
    locate("_", reverse(@replicate))-1 ) t0 ,
@t1 := LEFT(@replicate, length(@replicate) - locate("_", reverse(@replicate)) ) t1,
@t2 := LEFT(@t1, length(@t1) - locate("_", reverse(@t1)) ) t2, 
@t3 :=  locate("_", reverse(@replicate))  +  locate("_", reverse(@t1)) t3,
@t4 := CONCAT( 'Yeast' , RIGHT(@replicate,  @t3) ) as t4
#tina.sequence as sequence
#l.mrm_key as present_in_SRMAltas
#q3, Intensity, Ion_description, parent_charge, m.sequence as modseq, collision_energy,
#m.protein_name as mrmatlas_protein, tina.protein as protein
from hroest.tmp_tinasproteins tina
inner join hroest.tmp_tppsearch_analysis mann on mann.sequence = tina.sequence
order by sequence
;


select #count(*) occurence, 
protein, sds, fraction, replicate, fmap.id, fmap.file 
from ttmp tina
left join hroest.testfm fmap on fmap.file = t4
inner join hroest.manndata_lfq_sorted manndata on fmap.id = manndata.file
#inner join ddb.peptideOrganism po on po.peptide_key = tina.peptide_key
#left join hroest.MRMPepLink_final l on l.peptide_key = tina.peptide_key
#order by tina.protein, fraction, sds
#group by tina.protein, replicate
#
#group by tina.protein, replicate
#having occurence > 25
order by tina.protein#, occurence DESC
limit 20 ;


drop table hroest.testfm ;
create table hroest.testfm as 

select id, @replicate := file oldfile, 
@t1 := LEFT(@replicate, length(@replicate) - locate("_", reverse(@replicate)) ) t1,
@t2 := LEFT(@t1, length(@t1) - locate("_", reverse(@t1)) ) t2, 
@t3 :=  locate("_", reverse(@replicate))  +  locate("_", reverse(@t1)) t3,
@t4 := CONCAT( 'Yeast' , RIGHT(@replicate,  @t3) )  file
from hroest.manndata_lfq_filemapping 
;

select * from hroest.tmp_tppsearch_analysis_2 limit 2;


select count(*), sequence from 
hroest.tmp_tppsearch_analysis 
group by common_name
limit 2;






    }}}


}}}




KPADLASLLLNSAGDAQGDEAPALK                   
| YEL031W | LANVSAVTNIIR                    
| YEL031W | FQFSSALK                        
| YEL031W | GTVITPEIR                       
| YEL031W | SDDNQLLFR                       



select average, stdev, obs, file from manndata_lfq 
where peptide = 'SDDNQLLFR'
and average != 0
and obs = 3
;


select average, stdev, obs, file from manndata_lfq 
where peptide = 'GTVITPEIR'
and average != 0
#and obs = 3
;


}}}
#################################################################


c.execute()

###########################################################################
### SRM without RT information {{{ ##NO rt information
#

import numpy
query = """
select q1, q3 
from %(peptable)s 
inner join %(tratable)s on group_id = transition_group
where 
q1 between 400 and 1200 
and q3 between 400 and 1200 
and isotope_nr = 0
""" % {'peptable' : par.peptide_table ,
       'tratable' : par.transition_table }
cursor.execute( query)
lines = cursor.fetchall()

mybins = [ [0 for i in range(0, 1200) ] for i in range(0, 1200 / 0.7+ 1 ) ]
for l in lines:
    q1 = int( numpy.floor( l[0]/ 0.7 ))
    q3 = int( numpy.floor( l[1] ))
    mybins[q1][q3] += 1


total = 0
nr_bins = 0
mymax = 0
single_spots = 0
mybest = (0,0)
for out in range(400, 1200/0.7+1):
    for inn in range(400, 1200):
        total += mybins[out][inn]
        if mybins[out][inn] > mymax:
            mymax = mybins[out][inn]
            mybest = (out, inn)
        nr_bins += 1
        if mybins[out][inn] == 1: single_spots += 1

avg = 1.0* total / nr_bins
unique_pct = single_spots * 1.0 / total

}}}


{{{def get_max_charge(sequence):

def get_max_charge(sequence):
    #Calculate the max. charge state feasible for the peptide
    maxCharge = 1 #Any peptide may be charged at the N-terminus
    for aa in sequence:
        if aa in ['K', 'R', 'D', 'E', 'H']:
             maxCharge += 1
    return maxCharge

}}}



#################################################################
#{{{ 2010.12.14 Create mixture of N14 and N15

if True:
    #1000 entries / 9 s with fast (executemany), 10.5 with index
    #1000 entries / 31 s with normal (execute), 34 with index
    transition_table = 'hroest.srmTransitions_yeast_N14N15'
    peptide_table = 'hroest.srmPeptides_yeast_N14N15'
    c.execute('truncate table ' + peptide_table)
    c.execute('truncate table ' + transition_table)
    mass_bins = [ []  for i in range(0, 10000) ]
    rt_bins = [ []  for i in range(-100, 500) ]
    tmp_c  = db.cursor()
    progressm = progress.ProgressMeter(total=len(rows), unit='peptides')
    transition_group = 0
    for row in rows:
        progressm.update(1)
        transition_group += 1 
        #if ii >= 1000: break #####FOR TESTING ONLY
        for mycharge in [2,3]:  #precursor charge 2 and 3
            peptide = collider.get_peptide_from_table(t, row)
            peptide.charge = mycharge
            if modify_cysteins: peptide.modify_cysteins()
            peptide.create_fragmentation_pattern(R)
            if insert_db:
                #insert peptide into db
                collider.insert_peptide_in_db(peptide, db, peptide_table,
                                              transition_group=transition_group)
        #we want to insert the fragments only once per peptide
        if insert_db:
            #insert fragment charge 1 and 2 into database
            peptide.mass_H = R.mass_H
            collider.fast_insert_in_db(peptide, db, 1, transition_table,
                                       transition_group=transition_group)
            collider.fast_insert_in_db(peptide, db, 2, transition_table,
                                       transition_group=transition_group)
    ####
    #
    R.recalculate_monisotopic_data_for_N15()
    mass_bins = [ []  for i in range(0, 10000) ]
    rt_bins = [ []  for i in range(-100, 500) ]
    tmp_c  = db.cursor()
    progressm = progress.ProgressMeter(total=len(rows), unit='peptides')
    #
    #
    #do NOT reset transition_group
    for row in rows:
        progressm.update(1)
        transition_group += 1 
        #if ii >= 1000: break #####FOR TESTING ONLY
        for mycharge in [2,3]:  #precursor charge 2 and 3
            peptide = collider.get_peptide_from_table(t, row)
            peptide.charge = mycharge
            if modify_cysteins: peptide.modify_cysteins()
            peptide.create_fragmentation_pattern(R)
            if insert_db:
                #insert peptide into db
                collider.insert_peptide_in_db(peptide, db, peptide_table,
                                              transition_group=transition_group)
        #we want to insert the fragments only once per peptide
        if insert_db:
            #insert fragment charge 1 and 2 into database
            peptide.mass_H = R.mass_H
            collider.fast_insert_in_db(peptide, db, 1, transition_table,
                                       transition_group=transition_group)
            collider.fast_insert_in_db(peptide, db, 2, transition_table,
                                       transition_group=transition_group)





#A2 create the additional parent ions for the isotope patterns
PATTERNS_UP_TO = 3 #up to how many isotopes should be considered
cursor = c
vals = "peptide_key, q1_charge, q1, modified_sequence, ssrcalc, isotope_nr, transition_group"
query = "SELECT parent_id, %s FROM %s" % (vals, peptide_table)
cursor.execute( query )
allpeptides =  cursor.fetchall()
prepared = []
for p in allpeptides:
    q1_charge = p[2]
    for i in range(1,1+PATTERNS_UP_TO):
        prepared.append( [p[1], p[2], p[3] + (R.mass_diffC13 * i* 1.0) / q1_charge,
                              p[4], p[5], i, p[7]] )

q = "INSERT INTO %s (%s)" % (peptide_table, vals)  + \
     " VALUES (" + "%s," *6 + "%s)"
c.executemany( q, prepared)




}}}
#################################################################




Some SQL cleanup stuff, find big tables
\url{http://www.mysqlperformanceblog.com/2008/03/17/researching-your-mysql-table-sizes/}

select TABLE_NAME, TABLE_SCHEMA,
concat(round(data_length/(1024*1024*1024),2),'G') DATA,
concat(round(INDEX_LENGTH/(1024*1024*1024),2),'G') IDX
from information_schema.TABLES 
where TABLE_SCHEMA = 'hroest'
#order by data_length ;

select TABLE_NAME, TABLE_SCHEMA,
concat(round(data_length/(1024*1024),2),'M') DATA,
concat(round(INDEX_LENGTH/(1024*1024),2),'M') IDX
from information_schema.TABLES 
where TABLE_SCHEMA = 'hroest'
#order by data_length ;


SELECT count(*) TABLES,
concat(round(sum(table_rows)/1000000,2),'M') rows,
concat(round(sum(data_length)/(1024*1024*1024),2),'G') DATA,
concat(round(sum(index_length)/(1024*1024*1024),2),'G') idx,
concat(round(sum(data_length+index_length)/(1024*1024*1024),2),'G') total_size,
round(sum(index_length)/sum(data_length),2) idxfrac
FROM information_schema.TABLES
where TABLE_SCHEMA = 'hroest';





































truncate srmTransitions_perf_test;
insert into srmTransitions_perf_test (srm_id, group_id, q3_charge,
q3, type, fragment_number, q3_floor)
select 
srm_id, group_id, q3_charge,
q3, type, fragment_number,
floor(q3) 
from srmTransitions_yeast 
;



truncate srmPeptides_perf_test;
insert into srmPeptides_perf_test (parent_id, peptide_key, modified_sequence,
q1, ssrcalc, isotope_nr, transition_group, q1_charge, q1_floor, ssrcalc_floor) 
select parent_id, peptide_key,
modified_sequence, q1, ssrcalc, isotope_nr, transition_group, q1_charge,
floor(q1), floor(ssrcalc)
from srmPeptides_yeast 
;

select * from  srmPeptides_perf_test;

































 #why are they different?
 MRMAtlas_yeast_201009_final_top10                           |
 | MRMAtlas_yeast_201009_final_top25                           |
 | MRMAtlas_yeast_201009_final_top5                            |
 | MRMAtlas_yeast_qqq_201002_q2    



#has ions with more than one description / possibility
#unknown ions
#isotopes
#parent ions
#ions with square brackets [b9/-0.15]  => they appeared already?
select @seq := 'AAAAEKNVPLYKHLADLSK/3' ;
select @seq := 'DHITNAWHVPVTAQITEK/2' ;
select @seq := 'SSNSLDNQESSQQR/3 ' ;
select sequence, parent_mass, q3, Intensity, Ion_description
from MRMAtlas_yeast_qqq_201002_q2 where sequence =  @seq ;

select modified_sequence, q1, q3, Intensity, Ion_description
from MRMAtlas_yeast_201009_final_top25 where 
modified_sequence = @seq ;



select * from 
MRMAtlas_yeast_qqq_201002_q2 order by rand() limit 2;


def print_trans_collisions(par, db, p_id = 1, q3_low = 300, q3_high = 2000,
    query1_add = 'and q1_charge = 2 and q3_charge = 1'):
    res_str = ''
    cursor = db.cursor()
    non_unique = {}
    non_unique_q1 = {}
    non_unique_ppm = {}
    non_unique_clash = {}
    query1 = """
    select q3, srm_id, q1, ssrcalc, sequence, type, %(pep)s.peptide_key
    from %(pep)s
    inner join %(trans)s
      on %(pep)s.transition_group = %(trans)s.group_id
    inner join ddb.peptide on %(pep)s.peptide_key = peptide.id
    where parent_id = %(parent_id)s
    and q3 > %(q3_low)s and q3 < %(q3_high)s         %(query_add)s
    """ % { 'parent_id' : p_id, 'q3_low' : q3_low,
           'q3_high' : q3_high, 'query_add' : par.query1_add,
           'pep' : par.peptide_table, 'trans' : par.transition_table }
    nr_transitions = cursor.execute( query1 )
    transitions = cursor.fetchall()
    q1 = transitions[0][2]
    ssrcalc = transitions[0][3]
    sequence = transitions[0][4]
    peptide_key = transitions[0][6]
    prt = '\nAnalysing 2+ fragment ions of peptide nr %s\nWith sequence %s, q1 at %s and ssrcalc %s\n'
    res_str += prt % (p_id, sequence, q1, ssrcalc) 
    res_str += 'q3\tion      q3\tion q1   SSRcal dq3   sequence   q3charge\n'
    #
    query2 = """
    select q3, q1, ssrcalc, sequence, type, q3_charge, fragment_number
    from %(pep)s
    inner join %(trans)s
      on %(pep)s.transition_group = %(trans)s.group_id
    inner join ddb.peptide on %(pep)s.peptide_key = peptide.id
    where ssrcalc > %(ssrcalc)s - %(ssr_window)s 
        and ssrcalc < %(ssrcalc)s + %(ssr_window)s
    and q1 > %(q1)s - %(q1_window)s and q1 < %(q1)s + %(q1_window)s
    and %(pep)s.peptide_key != %(peptide_key)d
    and q3 > %(q3_low)s and q3 < %(q3_high)s
    %(query_add)s
    """ % { 'q1' : q1, 'ssrcalc' : ssrcalc, 'peptide_key' : peptide_key,
           'q3_low':q3_low,'q3_high':q3_high, 'q1_window' : par.q1_window,
           'query_add' : par.query2_add, 'ssr_window' : par.ssrcalc_window,
           'pep' : par.peptide_table, 'trans' : par.transition_table }
    tmp = cursor.execute( query2 )
    collisions = cursor.fetchall()
    #here we loop through all possible combinations of transitions and
    #potential collisions and check whether we really have a collision
    q3_window_used = par.q3_window
    isy = True; jj = len( transitions ) / 2 + 1
    jj = 0
    for ii,t in enumerate(transitions):
        ##before we didnt have y and b ions labelled
        if ii == len( transitions ) / 2: res_str += ( "=" * 75 + '\n'); isy = False; jj = 0
        #if isy: jj -= 1
        #else:  jj += 1
        jj += 1
        if par.ppm: q3_window_used = par.q3_window * 10**(-6) * t[0]
        this_min = q3_window_used + 1
        for c in collisions:
            if abs( t[0] - c[0] ) <= this_min:
                this_min = abs( t[0] - c[0] )
                non_unique[ t[1] ] = c
                non_unique_q1[ t[1] ] = q1 - c[1]
                non_unique_ppm[ t[1] ] = (t[0] - c[0] ) * 10**6 / t[0]
        if this_min < q3_window_used:
            c = non_unique[ t[1] ]
            res_str += ("""%.2f\t%s%s <==> %.2f\t%s%s %07.2f %.2f %.2f %s %s\n""" %
                  (t[0], t[5], jj, c[0], c[4], c[6], c[1], c[2], t[0] - c[0], c[3],
                   c[5] ) )
        else:
            res_str += ("""%.2f\t%s%s \n""" % (t[0], t[5], jj) )
    return res_str


class TransitionPeptide(object):
    def __init__(self):
        self.transitions = []
    def someprint(self):
        res_str = ''
        prt = '\nAnalysing peptide nr %s\nWith sequence %s, q1 at %s and ssrcalc %s\n'
        res_str += prt % (self.p_id, self.sequence, self.q1, self.ssrcalc) 
        res_str += 'q3\tion      q3\tion q1   SSRcal dq3   sequence   q3charge\n'
        return res_str
    def latex_describe(self):
        prt = '\nAnalysing peptide nr %s\nWith sequence %s, q1 at %s and ssrcalc %s\n'
        prt = 'Analysing peptide %s \\\\ With sequence %s, q1 at %s and ssrcalc %s'
        return prt % (self.p_id, self.sequence, self.q1, self.ssrcalc) 
