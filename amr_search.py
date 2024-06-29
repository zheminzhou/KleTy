from kleborate.resBLAST import resblast_one_assembly, read_class_file
try :
    from configure import dbs
except :
    from .configure import dbs


dbs = dbs['amr_dbs']

drug_types = [  ["Bla_acquired","AMR:BETA-LACTAM"],
                ["Col_acquired","AMR:COLISTIN"],
                ["Flq_acquired","AMR:QUINOLONE"],
                ["Bla_ESBL_inhR_acquired","AMR:INHIBITOR-RESISTANT"],
                ["Bla_inhR_acquired","AMR:INHIBITOR-RESISTANT"],
                ["AGly_acquired","AMR:AMINOGLYCOSIDE"],
                ["Bla_Carb_acquired","AMR:CARBAPENEM"],
                ["Bla_ESBL_acquired","AMR:ESBL"],
                ["Fcyn_acquired","AMR:FOSFOMYCIN"],
                ["MLS_acquired","AMR:MACROLIDE"],
                ["Phe_acquired","AMR:PHENICOL"],
                ["Rif_acquired","AMR:RIFAMYCIN"],
                ["Gly_acquired","AMR:GLYCOPEPTIDES"],
                ["Sul_acquired","AMR:SULFONAMIDE"],
                ["Tet_acquired","AMR:TETRACYCLINE"],
                ["Tgc_acquired","AMR:TIGECYCLINE"],
                ["Tmt_acquired","AMR:TRIMETHOPRIM"],
                ["Col_mutations","AMR:COLISTIN"],
                ["Flq_mutations","AMR:QUINOLONE"],
                ["Omp_mutations","AMR:BETA-LACTAM"],
                ["Bla_chr","AMR:BLA_INTRINSIC"], 
                # ["SHV_mutations","AMR:BLA_INTRINSIC"],
            ]



def amr_search(contigs, n_proc) :
    gene_info, _, _ = read_class_file(dbs['amr_class'])
    hits = resblast_one_assembly(contigs, gene_info, dbs['qrdr'], dbs['trunc'], dbs['omp'], dbs['CARD'], 60., 90., 60., 90.)
    
    # AMR:AMINOGLYCOSIDE AGly_acquired
    # AMR:BETA-LACTAM Bla_acquired Bla_chr Omp_mutations
    # AMR:CARBAPENEM Bla_Carb_acquired
    # AMR:ESBL Bla_ESBL_acquired
    # AMR:INHIBITOR-RESISTANT Bla_ESBL_inhR_acquired SHV_mutations
    # AMR:COLISTIN  Col_acquired Col_mutations
    # AMR:FOSFOMYCIN Fcyn_acquired
    # AMR:MACROLIDE MLS_acquired
    # AMR:PHENICOL Phe_acquired
    # AMR:QUINOLONE Flq_acquired Flq_mutations
    # AMR:RIFAMYCIN Rif_acquired
    # AMR:GLYCOPEPTIDES Gly_acquired
    # AMR:SULFONAMIDE Sul_acquired
    # AMR:TETRACYCLINE Tet_acquired
    # AMR:TIGECYCLINE Tgc_acquired
    # AMR:TRIMETHOPRIM Tmt_acquired
    # AMR:BLA_INTRINSIC Bla_chr

    # AMR:BLEOMYCIN|STRESS:COPPER|STRESS:MERCURY|STRESS:NICKEL|STRESS:SILVER|STRESS:TELLURIUM|STRESS:ARSENIC|STRESS:FLUORIDE|STRESS:QUATERNARY_AMMONIUM|VIRULENCE:clb|VIRULENCE:iro|VIRULENCE:iuc|VIRULENCE:rmp|VIRULENCE:ybt|Others|REPLICON:INC_TYPE|REPLICON:MOB_TYPE|REPLICON:MPF_TYPE

    res = {d:[] for g, d in drug_types }
    for g, d in drug_types :
        res[d].extend(hits.get(g, []))
    return res
