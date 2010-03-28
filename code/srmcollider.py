import MySQLdb
import sys 
sys.path.append( '/home/hroest/msa/code' )
from utils_h import utils
db = MySQLdb.connect(read_default_file="~/.my.cnf")
c = db.cursor()
c2 = db.cursor()

def get_collisions( q1, q3, ssrcalc, q1_window, q3_window, ssrcalc_window, 
                   exp_key):
    select_temp = """
    select * from hroest.srmClashes 
    inner join ddb.peptide on peptide.id = srmClashes.peptide_key
    where q1 between %s and %s
    and q3 between %s and %s
    and ssrcalc between %s and %s
    and srmClashes.experiment_key = %s
    ;
    """
    q1_low      =  q1 - q1_window /2.0
    q1_high     =  q1 + q1_window /2.0
    q3_low      =  q3 - q3_window /2.0   
    q3_high     =  q3 + q3_window /2.0
    ssrcalc_low =  ssrcalc - ssrcalc_window /2.0 
    ssrcalc_high=  ssrcalc + ssrcalc_window /2.0   
    ##q1_low      =   518
    ##q1_high     =   519
    ##q3_low      =   740
    ##q3_high     =   760
    ##ssrcalc_low =   20
    ##ssrcalc_high=   25
    t = utils.db_table( c2 )
    t.read( select_temp, [q1_low, q1_high, q3_low, q3_high, 
                            ssrcalc_low, ssrcalc_high, exp_key])
    #TODO implement
    #res = t.fetchall()
    res = []
    while True:
        r = t.fetchone()
        if r == None: break
        res.append( r)
    return t, res

input = """
DDGSGVDIIDRPSMCLEYTTSK   
DLEILPAGDLTEIGEK         
AVGIGFIAVGIIGYAIK        
DQTSNLR                  
DQLIENCSK                
NFMPLEAK                 
DLTETCR                  
CNATPETFLQDIDK           
SLGQELIR                 
DGISQEGTIENAPANPSK       
SFLVILPAGLK              
TQMSHTLCVTTQLPR          
DIFTTLVSPELLSTCK         
GILQYVWFKPFYCFGTLICSAWK  
ENGFSMFK                 
HLQGPSDLVK               
QTHPVIIPSFASWFDISK       
IEDDLLFR                 
SGLVPAQFIEPVR            
AISSKPLSVR               
"""
q1_w = 1
q3_w = 1
ssr_w = 4
exp_key = 3120 #yeast

#sanitize input
import re
seqs = "'"
for inp in input.split():
    #only alphanumeric
    seqs += filter(str.isalnum, inp)+ "','"

seqs = seqs[:-2]

#get input transitions
c.execute( 'use ddb;')
c2.execute( 'use ddb;')
input_q  = """
select * from peptide
inner join peptideOrganism on peptideOrganism.peptide_key = peptide.id
inner join hroest.srmClashes on srmClashes.peptide_key = peptide.id
where peptide.experiment_key = %s
and peptide.sequence in 
(%s)
""" % (exp_key, seqs)
t = utils.db_table( c2 )
t.read( input_q )
result = t.fetchall_groupBy('sequence')

for sequence in result.values():
    useable = []
    collisions= []
    for row in sequence:
        q1  = t.row( row, 'q1')
        q3  = t.row( row, 'q3')
        ssr = t.row( row, 'ssrcalc')
        occurence = t.row( row, 'genome_occurence')
        mysequence = t.row( row, 'sequence')
        tt, rr = get_collisions( q1, q3, ssr, q1_w, q3_w, ssr_w, exp_key)
        assert len( rr) > 0  #this would mean the sequence is not in the DB
        if len( rr) == 1: useable.append( row )
        else: 
            collisions.append( [row, [rrr for rrr in rr 
                  if tt.row(rrr, 'sequence') != mysequence]  ] )

    print mysequence, q1
    print "useable:", [ (t.row(u, 'type'), t.row(u, 'q3')) for u in useable]
    print "collisions: \n"
    for col in collisions:
        myrow = col[0]
        print t.row( myrow, 'type'), '::::', [( tt.row(c, 'sequence') , 
               tt.row(c, 'type'), 
               tt.row(c, 'q3'), tt.row(c, 'ssrcalc') ) for c in col[1] ] 


###########################################################################
###########################################################################
#                              SRM Collider                               #
###########################################################################
###########################################################################

