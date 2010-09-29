#to correleate ssrcalc with RT
select distinct
Tr_original, ssrcalc
from hroest.srmPeptides_yeast_allCAM
inner join hroest.srmTransitions_yeast_allCAM
ON parent_key = parent_id
inner join hroest.yeast_dp_data_light ON 
yeast_dp_data_light.modified_sequence = srmPeptides_yeast_allCAM.modified_sequence
where 
    prec_z = q1_charge
and frg_type = type
and frg_z = q3_charge
and frg_nr = fragment_number
;
