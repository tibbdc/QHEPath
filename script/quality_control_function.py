import cobra
import re
import math
import pandas as pd
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import pfba
from cobra.io.dict import model_to_dict, model_from_dict, metabolite_from_dict, gene_from_dict, reaction_from_dict


def add_rxn(model,ri,rr):
    """Adding reactions to the model.

    Args:
    
    * model (cobra.Model): A Model object
    * ri (str): ID of reaction
    * rr (str): ID of reaction equation

    """
    reaction = Reaction(ri)
    model.add_reactions([reaction])
    reaction.build_reaction_from_string(rr)

def rxn_no_annotation(model,rid):
    """Determine whether there are annotations for a reaction in the model.

    Args:
    
    * model (cobra.Model): A Model object
    * rid (str): ID of reaction
    
    return: If there is no annotation information in the reaction, return True.

    """
    rxn = model.reactions.get_by_id(rid)
    if rxn.annotation == []:
        return True
    if len(rxn.annotation) == 1 and  'MetaNetX (MNX) Equation' in str(rxn.annotation): 
        return True

def find_infinite_rxn(model,objective_reaction,constraint_value,rxn_grate):
    """Iteratively exclude the reaction associated with the highest penalty in each pFBA computation result 
        to eliminate specific types of errors within the model.

    Args:
    
    * model (cobra.Model): A Model object
    * objective_reaction(str): ID of objective reaction
    * constraint_value(float): The threshold for different error
    * rxn_grate (dict): Reactions in the model and their corresponding penalties
    
    return: List of closed reactions.

    """ 
    grate_delete_rxn = []
    n=0
    while (model.optimize().fluxes[objective_reaction] - constraint_value)>1e-6:
        n=n+1
        print(n)
        model.objective = objective_reaction
        pfba_solution = pfba(model)
        need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
        fba_rxn_dic={}
        for i in need_fluxes.index:
                fba_rxn_dic[i] = rxn_grate[i]
        for key,value in fba_rxn_dic.items():
            if value == max(fba_rxn_dic.values()): # 设置罚分最高的反应的bian'j
                if fba_rxn_dic[key]!= 0:
                    model.reactions.get_by_id(key).bounds=(0,0)
                    grate_delete_rxn.append(key)
                    print(model.reactions.get_by_id(key),fba_rxn_dic[key])
    return grate_delete_rxn

# 反应加回去，定位出造成无限生成的反应
def locating_rxn(model,check_model,objective_reaction,constraint_value,grate_delete_rxn,outputfile,unbalance_rxn):
    """Reintroduce reactions with high penalties back into the universal model to locate the reactions truly causing errors.

    Args:
    
    * model (cobra.Model): A Model object,the universal model
    * check_model (cobra.Model): A Model object, the universal model after determining the directions
    * constraint_value(float): The threshold for different error
    * grate_delete_rxn(list): List of closed reactions
    * outputfile (txt): Pathway output in a text (txt) file
    * unbalance_rxn (list): List of reaction IDs with mass imbalance
    
    return: List of reactions causing errors.

    """ 
    wuxian_rxn=[]
    n=0
    for i in grate_delete_rxn:
        n= n+1
        print(i,n)
        model.reactions.get_by_id(i).bounds=check_model.reactions.get_by_id(i).bounds
        pfba_solution=pfba(model)
        if (pfba_solution.fluxes[objective_reaction]-constraint_value)>1e-6:
            need_fluxes =pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
            file= outputfile + i +'.txt'
            pfba_pathway_output(model,file,need_fluxes)
            
            rxn = model.reactions.get_by_id(i)
            if i in unbalance_rxn:
                rxn.bounds = (0,0)
                wuxian_rxn.append(i)
                print(rxn,'unbalance',rxn.check_mass_balance())
            else:
                #if rxn.bounds == (-1000,1000): # 反应为可逆
                if pfba_solution.fluxes[i] > 0:
                    model.reactions.get_by_id(i).bounds=(-1000,0)
                    wuxian_rxn.append(i)
                    print(rxn,'direction_forward',rxn.check_mass_balance())
                if pfba_solution.fluxes[i] < 0:
                    model.reactions.get_by_id(i).bounds=(0,1000)
                    wuxian_rxn.append(i)
                    print(rxn,'direction_backward',rxn.check_mass_balance())
    return wuxian_rxn 

def reduced_degree(model,metid):
    """Calculate the reduction degree of metabolites based on their formulas.

    Args:
    
    * model (cobra.Model): A Model object,the universal model
    * metid (str): ID of reaction
    
    return: the reduction degree of metabolites.

    """ 
    met = model.metabolites.get_by_id(metid)
    degree = 0
    for e,en in met.elements.items():
        if e == 'C':
            degree = degree  + 4*en
        if e == 'H':
            degree = degree  + 1*en
        if e == 'O':
            degree = degree  + (-2)*en
        if e == 'N':
            degree = degree  + (-3)*en
        if e == 'P':
            degree = degree  + (5)*en
        if e == 'S':
            degree = degree  + (6)*en
    degree = degree - met.charge
    return degree
    

def max_reduced_degree_rate(model,metid):
    """Calculate the maximum theoretical yield of reduction for metabolites when 10 mmol of glucose is used as the substrate.

    Args:
    
    * metid (str): ID of reaction
    
    return: the maximum theoretical yield of reduction for metabolites.

    """ 
    if round(reduced_degree(model,metid),5) != 0:

        max_rate = 10*(reduced_degree(model,'glc__D_e')/reduced_degree(model,metid))
    else:
        max_rate = 0
    return max_rate

def pfba_pathway_output(model,outputfile,fluxes):
    """Exporting pathways calculated by pFBA to a text (txt) file.

    Args:
    
    * model (cobra.Model): A Model object of universal model (CSMN)
    * outputfile(str): Output file path
    * fluxes(dataframe): DataFrame returned by pFBA calculation

    """
    flux = open(outputfile, 'w')
    for r,v in fluxes.iteritems():
        if abs(v)>1e-6:
            v=round(v,5)
            check = model.reactions.get_by_id(r).check_mass_balance()
            flux.write(r+'\t'+str(v)+'\t'+model.reactions.get_by_id(r).build_reaction_string()+'\t'+str(check)+'\n')
    flux.close()



def transport_rxn(model,rid):
    """Determining whether a reaction is an transport reaction.

    Args:
    
    * model (cobra.Model): A Model object of universal model
    * rid (str): ID of reaction
    
    return: if the reaction is a transport reaction return True
    
    """
    rxn = model.reactions.get_by_id(rid)
    reactants_mets = [str(m).split('_')[0] for m in rxn.reactants]
    products_mets = [str(m).split('_')[0] for m in rxn.products]
    if sorted(reactants_mets) == sorted(products_mets):
        return True
    

def exchange_rxn(model,rid):
    """Determining whether a reaction is an exchange reaction.

    Args:
    
    * model (cobra.Model): A Model object of universal model
    * rid (str): ID of reaction
    
    return: if the reaction is a exchange reaction return True
    
    """
    rxn = model.reactions.get_by_id(rid)
    products_mets = [str(m).split('_')[0] for m in rxn.products]
    if len(products_mets) == 0:
        return True


def unbalance_coefficient(model,unri,balance_rxn): 
    """Compare whether the metabolites are the same in balanced and unbalanced reactions.

    Args:
    
    * model (cobra.Model): A Model object of universal model
    * unri (str): ID of balanced reaction
    * balance_rxn (list): List of reaction IDs for balanced reactions
    
    return: Balanced reactions with metabolites identical to unbalanced reactions.
    
    """
    unr = model.reactions.get_by_id(unri)
    unm = [m.id for m in unr.metabolites]
    for balri in balance_rxn:
        if not transport_rxn(model,balri):
            balr = model.reactions.get_by_id(balri)
            balm = [m.id for m in balr.metabolites]
            if sorted(unm) == sorted(balm):
                return balr 
            
def unbalance_miss_cofactor(model,unri,balance_rxn): 
    """Determine if unbalanced reactions are missing cofactors.

    Args:
    
    * model (cobra.Model): A Model object of universal model
    * unri (str): ID of balanced reaction
    * balance_rxn (list): List of reaction IDs for balanced reactions
    
    return: Balanced reactions with additional cofactors compared to unbalanced reactions.
    
    """
    cofactor = ['h','h2o','nad','nadh','nadp','nadph','fad','fadh2']
    unr = model.reactions.get_by_id(unri)
    unbal_met = [m.id for m in unr.metabolites]
    for balri in balance_rxn:
        if not transport_rxn(model,balri):
            balr = model.reactions.get_by_id(balri)
            bal_met = [m.id for m in balr.metabolites]
            diff1 = set(unbal_met) - set(bal_met)
            diff2 = set(bal_met) - set(unbal_met)
            diff3 = [m.split('_')[0]  for m in list(diff1 | diff2)]
            if len(diff3)!= 0 and set(diff3).issubset(cofactor):
                return balr
            
def unbalance_error_met(model,unri,balance_rxn): 
    """Determine if unbalanced reactions contain incorrect metabolites.

    Args:
    
    * model (cobra.Model): A Model object of universal model
    * unri (str): ID of balanced reaction
    * balance_rxn (list): List of reaction IDs for balanced reactions
    
    return: Balanced reactions with correct metabolites compared to unbalanced reactions.
    
    """
    currency_met = ['h','h2o']
    currency_met2 = ['nad','nadh','nadp','nadph','pqb6s','pq']
    unr = model.reactions.get_by_id(unri)
    unbal_rxn_met1 = [m.id.split('_')[0] for m in unr.metabolites if m.id.split('_')[0] not in currency_met]
    unbal_reactant_met = [m.id for m in unr.reactants if m.id.split('_')[0] not in currency_met]
    unbal_product_met = [m.id for m in unr.products if m.id.split('_')[0] not in currency_met]
    for balri in balance_rxn:
        if not transport_rxn(model,balri) and not set(unbal_rxn_met1).issubset(currency_met2):
            balr = model.reactions.get_by_id(balri)
            bal_reactant_met = [m.id for m in balr.reactants if m.id.split('_')[0] not in currency_met]
            bal_product_met = [m.id for m in balr.products if m.id.split('_')[0] not in currency_met]

            if len((set(unbal_reactant_met) -set(bal_reactant_met))) == 1 and len((set(bal_reactant_met) -set(unbal_reactant_met))) and len((set(unbal_product_met) -set(bal_product_met))) == 0 and len((set(bal_product_met) -set(unbal_product_met))) == 0:
                return balr
            if len((set(unbal_product_met) -set(bal_product_met))) == 1 and len((set(unbal_product_met) -set(bal_product_met))) == 1 and len((set(unbal_reactant_met) -set(bal_reactant_met))) == 0 and len((set(bal_reactant_met) -set(unbal_reactant_met))) == 0 :
                return balr

            if len((set(unbal_reactant_met) -set(bal_reactant_met))) == 0 and len((set(bal_reactant_met) -set(unbal_reactant_met))) == 1 and len((set(unbal_product_met) -set(bal_product_met))) == 0 and len((set(bal_product_met) -set(unbal_product_met))) == 0:
                return balr
            if len((set(bal_product_met) -set(unbal_product_met))) == 0 and len((set(unbal_product_met) -set(bal_product_met))) == 1 and len((set(unbal_reactant_met) -set(bal_reactant_met)))  == 0 and len((set(bal_reactant_met) -set(unbal_reactant_met))) == 0:
                return balr
            
def o2_rule(model,detaG_correct):
    """Determine the direction of metabolites containing oxygen.
        For reactions containing oxygen, their directions were assumed to be that of the oxygen consumption.

    Args:
    
    * model (cobra.Model): A Model object of universal model
    * detaG_correct (list): List of reactions corrected based on thermodynamics
    
    return: List of reactions with corrected directions.
    
    """

    # Determine inorganic reactions where metabolites contain oxygen.
    O2_inorganic_rxn=[]
    o2_mets=['h2o2','o2','h','h2o']
    for rea in model.reactions:
        if not transport_rxn(model,rea.id) and not exchange_rxn(model,rea.id):
                n=0
                s=0
                for met in rea.metabolites:
                    s=s+1
                    if str(met).split('_')[0] in o2_mets:
                        n=n+1
                if s==n:
                    O2_inorganic_rxn.append(rea.id)

    # Reactions producing oxygen through superoxides and peroxides.
    # SPODMpp: 2.0 h_p + 2.0 o2s_p --> h2o2_p + o2_p
    # SOD_m: 2.0 sox_m --> h2o2_m + o2_m
    # GTHP_CAT: 2.0 gthrd_c + 3.0 h2o2_c --> gthox_c + 4.0 h2o_c + o2_c 
    superoxide_rxn = []
    for r in model.reactions:
        reactants_mets=[str(m).split('_')[0] for m in r.reactants]
        products_mets=[str(m).split('_')[0]  for m in r.products]
        if 'o2s' in reactants_mets or 'sox' in reactants_mets:
            if 'h2o2' in products_mets and 'o2' in products_mets:
                superoxide_rxn.append(r.id)
                
        if 'o2s' in products_mets or 'sox' in products_mets:
            if 'h2o2' in reactants_mets and 'o2' in reactants_mets:
                superoxide_rxn.append(r.id)
                
        if 'gthrd' in reactants_mets or 'h2o2' in reactants_mets:
            if 'gthox' in products_mets and 'o2' in products_mets:
                superoxide_rxn.append(r.id)
                r.bounds = (0,1000)
        
        if 'gthrd' in products_mets or 'h2o2' in products_mets:
            if 'gthox' in reactants_mets and 'o2' in reactants_mets:
                superoxide_rxn.append(r.id)
                r.bounds = (-1000,0)

    # Reactions generating oxygen through photosynthesis.
    photosynthetic_rxn = []
    for r in model.reactions:
        if re.search(r'Photosystem',r.name):
            photosynthetic_rxn.append(r.id)

    O2_rxn = []
    for r in model.reactions:
        if not transport_rxn(model,r.id) and not exchange_rxn(model,r.id) and r.id not in O2_inorganic_rxn and r.id not in superoxide_rxn and r.id not in photosynthetic_rxn:
            if r.id != 'SOPSI' and r.id != 'SOPSII':
                if r.id not in detaG_correct:
                    reactants_mets=[str(m).split('_')[0] for m in r.reactants]
                    products_mets=[str(m).split('_')[0]  for m in r.products]
                    if 'o2' in reactants_mets:
                        if r.bounds == (-1000,1000) or r.bounds == (-1000,0):
                            O2_rxn.append(r.id)
                            r.bounds = (0,1000)
                            
                    if 'o2' in products_mets:
                        if r.bounds == (-1000,1000) or r.bounds == (0,1000):
                            O2_rxn.append(r.id)
                            r.bounds = (-1000,0)         
    return O2_rxn

def nh4_rule(model,detaG_correct,O2_rxn):
    """Determine the direction of metabolites containing NH3/NH4+.
        most reactions produced  NH3/NH4+, except when it reacted with ATP, 2-oxoglutarate, chorismate, 
        5-phospho-alpha-D-ribose-1-diphosphate, and UDP-N-acetyl-beta-L-fucosamine.

    Args:
    
    * model (cobra.Model): A Model object of universal model
    * detaG_correct (list): List of reactions corrected based on thermodynamics
    * O2_rxn (list): List of reactions corrected based on o2_rule
    
    return: List of reactions with corrected directions.
    
    """

    nh4_inorganic_rxn=[] 
    for r in model.reactions:
        if not transport_rxn(model,r.id) and not exchange_rxn(model,r.id):
            reactants_mets=[str(m).split('_')[0] for m in r.reactants]
            products_mets=[str(m).split('_')[0]  for m in r.products]
            if 'nh3' in reactants_mets and 'nh4' in products_mets:
                nh4_inorganic_rxn.append(r.id)
                
            if 'nh3' in products_mets and 'nh4' in reactants_mets:
                nh4_inorganic_rxn.append(r.id)

    nh3_nh4_rxn = []
    for r in model.reactions:
        if not transport_rxn(model,r.id) and not exchange_rxn(model,r.id) and r.id not in nh4_inorganic_rxn and r.id not in O2_rxn:
            if r.id not in detaG_correct:
                reactants_mets=[str(m).split('_')[0] for m in r.reactants]
                products_mets=[str(m).split('_')[0]  for m in r.products]
                if 'nh4' in reactants_mets or 'nh3' in reactants_mets:
                    if 'atp' not in reactants_mets and 'akg' not in reactants_mets and 'chor' not in reactants_mets and 'prpp' not in reactants_mets and 'udpacblfuc' not in reactants_mets:
                        if r.bounds == (-1000,1000) or r.bounds == (0,1000):
                            r.bounds = (-1000,0) 
                            nh3_nh4_rxn.append(r.id)
                            # print(r)
                if 'nh4' in products_mets or 'nh3' in products_mets:
                    if 'atp' not in products_mets and 'akg' not in products_mets and 'chor' not in products_mets and 'prpp' not in products_mets and 'udpacblfuc' not in products_mets:
                        if r.bounds == (-1000,1000) or r.bounds == (-1000,0):
                            nh3_nh4_rxn.append(r.id)
                            r.bounds = (0,1000)
                            # print(r)
    return nh3_nh4_rxn

def co2_rule(model,detaG_correct):
    """Determine the direction of metabolites containing CO2.
        Most reactions cannot proceed in the direction of CO2 fixation, 
        except for naturally known carbon fixation reactions and reactions utilizing CO2 and high-energy substrates 
        such as phosphoenolpyruvate (PEP) and ATP.

    Args:
    
    * model (cobra.Model): A Model object of universal model
    * detaG_correct (list): List of reactions corrected based on thermodynamics
    
    return: List of reactions with corrected directions.
    
    """
    CO2_rxn_all = []
    for rxn in model.reactions:
        if rxn.id not in detaG_correct:
            rxn_mets = [m.id.split('_')[0] for m in rxn.metabolites]
            rmets = [m.id.split('_')[0] for m in rxn.reactants]
            pmets = [m.id.split('_')[0] for m in rxn.products]
            if not transport_rxn(model,rxn.id) and not exchange_rxn(model,rxn.id):
                if rxn.bounds == (-1000,1000):
                    if 'co2' in rxn_mets:
                        CO2_rxn_all.append(rxn.id)
                        
                if rxn.bounds == (0,1000):
                    if 'co2' in rmets:
                        CO2_rxn_all.append(rxn.id)
                        
                if rxn.bounds == (-1000,0):
                    if 'co2' in pmets:
                        CO2_rxn_all.append(rxn.id)
                        
    # Inorganic reactions containing CO2
    co2_inorganic_rxn=[]
    co2_mets=['co2','hco3','h2co3','h','h2o']
    for ri in CO2_rxn_all:
        rxn = model.reactions.get_by_id(ri)
        rxn_mets = [str(m).split('_')[0] for m in rxn.metabolites]
        if set(rxn_mets).issubset(co2_mets):
            co2_inorganic_rxn.append(ri)

    co2_atp_pep_rxn=[]
    for ri in CO2_rxn_all:
        rxn = model.reactions.get_by_id(ri)
        reactants_mets=[str(m).split('_')[0]  for m in rxn.reactants]
        products_mets=[str(m).split('_')[0]  for m in rxn.products]
        if 'co2' in reactants_mets:
            if 'atp' in reactants_mets or 'pep' in reactants_mets:
                if 'atp' not in products_mets and 'itp' not in products_mets:
                    co2_atp_pep_rxn.append(ri)
                    
        if 'co2' in products_mets:
            if 'atp' in products_mets or 'pep' in products_mets:
                if 'atp' not in reactants_mets and 'itp' not in reactants_mets:
                    co2_atp_pep_rxn.append(ri)
                    
    # Natural carbon fixation reactions
    natural_co2_rxn = []
    #  RBPC: co2_c + h2o_c + rb15bp_c --> 2.0 3pg_c + 2.0 h_c
    for ri in CO2_rxn_all:
        if ri not in co2_inorganic_rxn:
            rxn = model.reactions.get_by_id(ri)
            rxn_mets = [m.id.split('_')[0] for m in rxn.metabolites]
            reactants_mets=[str(m).split('_')[0]  for m in rxn.reactants]
            products_mets=[str(m).split('_')[0]  for m in rxn.products]

            if 'co2' in reactants_mets and 'rb15bp' in reactants_mets and '3pg' in products_mets:
                natural_co2_rxn.append(ri)
                
            if 'co2' in products_mets and 'rb15bp' in products_mets and '3pg' in reactants_mets:
                natural_co2_rxn.append(ri)
                
                
            if 'co2' in reactants_mets and 'for' in products_mets: 
                natural_co2_rxn.append(ri)
                
            if 'co2' in products_mets and 'for' in reactants_mets:
                natural_co2_rxn.append(ri)
            
            if 'co2' in reactants_mets and 'for' in products_mets: 
                natural_co2_rxn.append(ri)
                
            if 'co2' in products_mets and 'for' in reactants_mets:
                natural_co2_rxn.append(ri)
            natural_co2_rxn.append('CODH_ACS')
            
            # OOR2r: akg_c + coa_c + fdxo_42_c <=> co2_c + fdxr_42_c + h_c + succoa_c
            # ICDHyr: icit_c + nadp_c <=> akg_c + co2_c + nadph_c
            if 'co2' in rxn_mets and 'fdxo' in rxn_mets and 'akg' in rxn_mets and 'fdxr' in rxn_mets and 'succoa' in rxn_mets: 
                natural_co2_rxn.append(ri)

            if 'co2' in rxn_mets and 'icit' in rxn_mets and 'akg' in rxn_mets and 'nadph' in rxn_mets and 'nadp' in rxn_mets: 
                natural_co2_rxn.append(ri)

            #POR5: coa_c + 2.0 flxso_c + pyr_c <=> accoa_c + co2_c + 2.0 flxr_c + h_c
            #POR: coa_c + fdxo_42_c + pyr_c <=> accoa_c + co2_c + fdxr_42_c + h_c
            if 'co2' in rxn_mets and 'flxso' in rxn_mets and 'pyr' in rxn_mets and 'flxr' in rxn_mets and 'accoa' in rxn_mets: 
                natural_co2_rxn.append(ri)
                
            if 'co2' in rxn_mets and 'fdxo' in rxn_mets and 'pyr' in rxn_mets and 'fdxr' in rxn_mets and 'accoa' in rxn_mets: 
                natural_co2_rxn.append(ri)
                
    co2_rule_rxn = []
    for ri in list(set(CO2_rxn_all) - set(co2_inorganic_rxn) - set(co2_atp_pep_rxn) -set(natural_co2_rxn)):
        if ri not in detaG_correct:
            rxn = model.reactions.get_by_id(ri)
            rxn_mets = [m.id.split('_')[0] for m in rxn.metabolites]
            reactants_mets=[str(m).split('_')[0]  for m in rxn.reactants]
            products_mets=[str(m).split('_')[0]  for m in rxn.products]
            if 'co2' in reactants_mets:
                if rxn.bounds ==(-1000,1000) or rxn.bounds ==(0,1000):
                    co2_rule_rxn.append(ri)
                    rxn.bounds =(-1000,0)
                    
            if 'co2' in products_mets:
                if rxn.bounds == (-1000,1000) or rxn.bounds == (-1000,0):
                    co2_rule_rxn.append(ri)
                    rxn.bounds = (0,1000)
                    
    return co2_rule_rxn

def atp_rule(model,detaG_correct):
    """Determine the direction of metabolites containing ATP, GTP, ITP, UTP.
        Most reactions involving ATP, GTP, ITP, UTP are in the consumption direction.

    Args:
    
    * model (cobra.Model): A Model object of universal model
    * detaG_correct (list): List of reactions corrected based on thermodynamics
    
    return: List of reactions with corrected directions.
    
    """
    # Proton-driven reactions
    ATP_synthase_rxn=[]
    ATP_list=['atp','adp','h','pi','h2o','h']
    for r in model.reactions:
        ATP_mets=[]
        for met in r.metabolites:
            if str(met)[-2:-1]=='_':
                ATP_mets.append(str(met)[:-2])
            elif str(met)[-3:-2]=='_':
                ATP_mets.append(str(met)[:-3])   
        if sorted(ATP_mets)==sorted(ATP_list):
            ATP_synthase_rxn.append(r.id)
    
    ATP_pi=[]
    ylcoa_met = ['accoa','succoa']    
    for r in model.reactions:
        if not transport_rxn(model,r.id) and not exchange_rxn(model,r.id) and r.id not in ATP_synthase_rxn:
            if r.id not in detaG_correct:
                reactants_mets=[str(m).split('_')[0] for m in r.reactants]
                products_mets=[str(m).split('_')[0]  for m in r.products]
                if 'atp' in reactants_mets or 'gtp' in reactants_mets or 'ctp' in reactants_mets or 'utp' in reactants_mets or 'itp' in reactants_mets:
                    if 'pi' in products_mets or 'ppi' in products_mets or 'pppi' in products_mets:
                        if len(set(products_mets) & set(ylcoa_met)) == 0 and r.lower_bound==-1000:  
                            r.bounds = (0,1000)
                            ATP_pi.append(r.id)
                            
                if 'atp' in products_mets or 'gtp' in products_mets or 'ctp' in products_mets or 'utp' in products_mets or 'itp' in products_mets:
                    if 'pi' in reactants_mets or 'ppi' in reactants_mets or 'pppi' in reactants_mets:
                        if len(set(reactants_mets) & set(ylcoa_met)) == 0 and r.lower_bound==-1000:   
                            ATP_pi.append(r.id)
                            r.bounds = (-1000,0)
    return ATP_pi

def ppi_h2o_rule(model): 
    """Determine the direction of metabolites containing triphosphate and diphosphate.
        The hydrolysis reaction of polyphosphoric acid is in the direction of hydrolysis.

    Args:
    
    * model (cobra.Model): A Model object of universal model
    
    return: List of reactions with corrected directions.
    
    """
    ppi_h2o=[] 
    for r in model.reactions:
        if not exchange_rxn(model,r.id):
            reactants_mets=[str(r_m).split('_')[0] for r_m in r.reactants]
            products_mets=[str(p_m).split('_')[0] for p_m in r.products]
            if 'h2o' in reactants_mets:
                for r_m in r.reactants:
                    if re.search(r'triphosphate|diphosphate',r_m.name) and r.lower_bound==-1000: 
                        if 'pi' in products_mets or 'ppi' in products_mets:
                            ppi_h2o.append(r.id)
                            r.bounds=(0,1000)
                            
            if 'h2o' in products_mets:            
                for p_m in r.products:
                    if re.search(r'triphosphate|diphosphate',p_m.name) and r.upper_bound==1000: 
                        if 'pi' in reactants_mets or 'ppi' in reactants_mets:
                            ppi_h2o.append(r.id)
                            r.bounds=(-1000,0)    
    return ppi_h2o

def ylcoa_h2o(model): 
    """Determine the direction of metabolites containing acyl-CoA.
        The hydrolysis of acyl-CoA to form CoA is the direction of hydrolysis.

    Args:
    
    * model (cobra.Model): A Model object of universal model
    * detaG_correct (list): List of reactions corrected based on thermodynamics
    
    return: List of reactions with corrected directions.
    
    """
    # Metabolites containing acyl-CoA
    ylcoa_met = []
    for met in model.metabolites:
        if re.search(r'yl-CoA|yl CoA|yl coenzyme A|yl Coenzyme A',met.name):
            ylcoa_met.append(str(met.id).split('_')[0])

    acyl_h2o = []
    for r in model.reactions:
        if not exchange_rxn(model,r.id):
            reactants_mets=[str(m).split('_')[0] for m in r.reactants]
            products_mets=[str(m).split('_')[0]  for m in r.products]
            if 'h2o' in reactants_mets and  len(set(ylcoa_met) & set(reactants_mets)) != 0 and len(reactants_mets) ==2 and 'coa' in products_mets:
                if r.bounds == (-1000,1000) or r.bounds == (-1000,0):
                    acyl_h2o.append(r.id)
                    r.bounds = (0,1000)
                                
            if 'h2o' in products_mets and len(set(ylcoa_met) & set(products_mets))!= 0 and len(products_mets) == 2 and 'coa' in reactants_mets:
                if r.bounds == (-1000,1000) or r.bounds == (0,1000):
                    acyl_h2o.append(r.id)
                    r.bounds = (-1000,0)
    return acyl_h2o

def ac_rule(model):
    """Determine the direction of metabolites containing high-energy compounds (ATP, Acyl-CoA).
        In addition to the reaction of acetate with high-energy compounds (ATP, Acyl-CoA), 
        the reaction containing acetate is the direction of producing acetate.

    Args:
    
    * model (cobra.Model): A Model object of universal model
    
    return: List of reactions with corrected directions.
    
    """

    ylcoa_met = []
    for met in model.metabolites:
        if re.search(r'yl-CoA|yl CoA|yl coenzyme A|yl Coenzyme A',met.name):
            ylcoa_met.append(str(met.id).split('_')[0])
            
    ac_rxn = []
    for r in model.reactions:
        if not transport_rxn(model,r.id) and not exchange_rxn(model,r.id):
            reactants_mets=[str(m).split('_')[0] for m in r.reactants]
            products_mets=[str(m).split('_')[0]  for m in r.products]
            if 'ac' in products_mets:
                if 'atp' not in products_mets and len(set(products_mets) & set(ylcoa_met)) == 0:
                    if r.bounds == (-1000,1000) or r.bounds == (-1000,0):
                        r.bounds = (0,1000)
                        ac_rxn.append(r.id)
                        # print(r)
            if 'ac' in reactants_mets:
                if 'atp' not in reactants_mets and len(set(reactants_mets) & set(ylcoa_met)) == 0:
                    if r.bounds == (-1000,1000) or r.bounds == (0,1000):
                        r.bounds = (-1000,0)
                        ac_rxn.append(r.id)
    return ac_rxn

def aldehyde_rule(model):
    """Determine the direction of metabolites containing aldehydes.
        The oxidation of aldehydes to acids by NA(D)P and ferredoxin is irreversible

    Args:
    
    * model (cobra.Model): A Model object of universal model
    
    return: List of reactions with corrected directions.
    
    """
    ate_met = []  
    for m in model.metabolites:
        if m.name.endswith('ate') or m.name.split(' (')[0].endswith('ate') or m.name.endswith('Acid') or  m.name.endswith('acid'):
            ate_met.append(m.id.split('_')[0])
    aldehyde_met = [] 
    for m in model.metabolites:
        if m.name.endswith('aldehyde') or m.name.endswith('anal') or re.search(r'Erythrose|erythrose',m.name):
            aldehyde_met.append(m.id.split('_')[0])

    aldehyde_ate_rxn = []        
    for r in model.reactions:
        if not exchange_rxn(model,r.id):
            reactants_mets=[str(m).split('_')[0] for m in r.reactants]
            products_mets=[str(m).split('_')[0]  for m in r.products]
            if 'nad' in reactants_mets or 'nadp' in reactants_mets or 'fdxo'in reactants_mets:
                if len(set(reactants_mets) & set(aldehyde_met)) != 0 and 'h2o' in reactants_mets and len(set(products_mets) & set(ate_met)) != 0 and 'o2' not in products_mets:
                    if r.bounds == (-1000,1000) or r.bounds == (-1000,0):
                        aldehyde_ate_rxn.append(r.id)
                        r.bounds = (0,1000)
                        
            if 'nad' in products_mets or 'nadp' in products_mets or 'fdxo'in products_mets:
                if len(set(products_mets) & set(aldehyde_met)) != 0 and 'h2o' in products_mets and len(set(reactants_mets) & set(ate_met)) != 0 and 'o2' not in reactants_mets:
                    if r.bounds == (-1000,1000) or r.bounds == (0,1000):
                        aldehyde_ate_rxn.append(r.id)
                        r.bounds = (-1000,0)
    return aldehyde_ate_rxn

def modify_direction_compartment(model):
    rxncl = []
    rxn_otherl = []

    for rxn in model.reactions:
        if not transport_rxn(model,rxn.id) and not exchange_rxn(model,rxn.id):
            rxn_mets = [m.id.split('_')[0] for m in rxn.metabolites]
            rxn_mets_c = [m.id for m in rxn.metabolites if m.id.split('_')[-1]=='c']
            if len(rxn_mets) == len(rxn_mets_c):
                rxncl.append(rxn.id)
            else:
                rxn_otherl.append(rxn.id)

    rxnc_direction =[]
    rxno_direction = []
    for rci in rxncl:
        rxnc = model.reactions.get_by_id(rci)
        rxnc_mets = [m.id.split('_')[0] for m in rxnc.metabolites]
        rmetc = [m.id.split('_')[0] for m in rxnc.reactants]
        pmetc = [m.id.split('_')[0] for m in rxnc.products]
        for roi in rxn_otherl:
            rxno = model.reactions.get_by_id(roi)
            rxno_mets = [m.id.split('_')[0] for m in rxno.metabolites]
            rmeto = [m.id.split('_')[0] for m in rxno.reactants]
            pmeto = [m.id.split('_')[0] for m in rxno.products]
            if set(rxnc_mets) == set(rxno_mets):
        
                if set(rmetc) == set(rmeto) and set(pmetc) == set(pmeto):
                    if rxnc.bounds != rxno.bounds:
                        rxno.bounds = rxnc.bounds
                        rxnc_direction.append(rxnc.id)
                        rxno_direction.append(rxno.id)
                            
                if set(rmetc) == set(pmeto) and set(pmetc) == set(rmeto):
                    if rxnc.bounds == (-1000,1000) and  rxno.bounds != (-1000,1000):
                        rxno.bounds = (-1000,1000)
                        rxnc_direction.append(rxnc.id)
                        rxno_direction.append(rxno.id)
                    if rxnc.bounds == (0,1000):
                        rxno.bounds = (-1000,0)
                        rxnc_direction.append(rxnc.id)
                        rxno_direction.append(rxno.id)
        
                    if rxnc.bounds == (-1000,0):
                        rxno.bounds = (0,1000)
                        rxnc_direction.append(rxnc.id)
                        rxno_direction.append(rxno.id)
                        
    return rxnc_direction

def check_rxn_balance_number(model):
    
    balance_rxn = []
    chargenan_rxn = []
    H_balance_rxn = []
    H2O_balance_rxn = []
    unbalance_rxn = []

    for rxn in model.reactions:
        if not exchange_rxn(model,rxn.id):
            rxn_balance_dic = rxn.check_mass_balance()
            if rxn_balance_dic == {}:
                balance_rxn.append(rxn.id)
            elif 'charge' in rxn_balance_dic.keys() and len(rxn_balance_dic) == 1 and math.isnan(rxn_balance_dic['charge']):  
                balance_rxn.append(rxn.id)
                chargenan_rxn.append(rxn.id) 
            elif "'H': " in str(rxn_balance_dic) and "'charge': " in str(rxn_balance_dic) and len(rxn_balance_dic) == 2 and rxn_balance_dic['charge'] == rxn_balance_dic['H'] and rxn_balance_dic['charge'] * rxn_balance_dic['H'] >=0: #{'charge': 1.0, 'H': 1.0}，电荷数和氢数相同，并且同号，认为此反应正确，不进行人为检查
                    H_balance_rxn.append(rxn.id)
                    balance_rxn.append(rxn.id)
            
            elif "'H': " in str(rxn_balance_dic) and "'O': " in str(rxn_balance_dic) and len(rxn_balance_dic) == 2 and abs(rxn_balance_dic['O']) == 1 and  abs(rxn_balance_dic['H']) == 2 and rxn_balance_dic['O'] * rxn_balance_dic['H'] >=0: #{'charge': 1.0, 'H': 1.0}，电荷数和氢数相同，并且同号，认为此反应正确，不进行人为检查
                H2O_balance_rxn.append(rxn.id)
                balance_rxn.append(rxn.id)
            else:
                unbalance_rxn.append(rxn.id)
    return balance_rxn,chargenan_rxn,H_balance_rxn,H2O_balance_rxn,unbalance_rxn

def check_carbon_production(model,element_id,rxn_grate,check_model,outputfile,unbalance_rxn):

    
    model.reactions.get_by_id('EX_'+element_id+'_e').bounds = (0,0)
    
    medium_rxn = ['EX_co2_e','EX_o2_e','EX_h_e','EX_h2o_e','EX_nh4_e','EX_pi_e','EX_so4_e','EX_fe3_e', 'EX_mn2_e','EX_fe2_e','EX_zn2_e', 'EX_mg2_e', 'EX_ca2_e', 'EX_ni2_e', 'EX_cu2_e', 'EX_cobalt2_e', 'EX_sel_e','EX_mobd_e', 'EX_k_e', 'EX_na1_e','EX_cl_e', 'EX_tungs_e', 'EX_slnt_e']
    for ri in medium_rxn:
        if ri!= 'EX_'+element_id+'_e':
            rxn = model.reactions.get_by_id(ri)
            rxn.bounds = (-1000,1000)
        
    met_wuxian = []
    with model:
        for met in model.metabolites:
            if re.search("C[A-Z]",str(met.formula)) or re.search("C[\d]",str(met.formula)): # 判断代谢物是否含碳，C后边是大写字母或者数字
                target_product = met.id
                objective_reaction = 'ADD_'+target_product
                if objective_reaction in model.reactions:
                    model.objective = objective_reaction
                    rxn_grate[objective_reaction] = 0
                else:
                    add_rxn(model,objective_reaction,target_product + ' --> ') #加入目标函数反应
                    model.objective = objective_reaction
                    rxn_grate[objective_reaction] = 0
                solution = model.optimize()
                print(target_product,solution.fluxes[objective_reaction])

                constraint_value = 1e-6
                if solution.fluxes[objective_reaction] > 1e-6:
                    met_wuxian.append(met.id)
                    grate_delete_rxn = find_infinite_rxn(model,objective_reaction,constraint_value,rxn_grate)
                    c_met = locating_rxn(model,check_model,objective_reaction,constraint_value,grate_delete_rxn,outputfile,unbalance_rxn)
                    c_met +=c_met
    return c_met

def check_other_elment_production(model,element_id,rxn_grate,check_model,outputfile,unbalance_rxn):
    model.reactions.get_by_id('EX_'+element_id+'_e').bounds = (0,0)
    
    medium_rxn = ['EX_co2_e','EX_o2_e','EX_h_e','EX_h2o_e','EX_nh4_e','EX_pi_e','EX_so4_e','EX_fe3_e', 'EX_mn2_e','EX_fe2_e','EX_zn2_e', 'EX_mg2_e', 'EX_ca2_e', 'EX_ni2_e', 'EX_cu2_e', 'EX_cobalt2_e', 'EX_sel_e','EX_mobd_e', 'EX_k_e', 'EX_na1_e','EX_cl_e', 'EX_tungs_e', 'EX_slnt_e']
    for ri in medium_rxn:
        if ri!= 'EX_'+element_id+'_e':
            rxn = model.reactions.get_by_id(ri)
            rxn.bounds = (-1000,1000)
    with model:    
        met_wuxian = []
        met = model.metabolites.get_by_id(element_id+'_c')
        target_product = met.id
        objective_reaction = 'ADD_'+target_product
        if objective_reaction in model.reactions:
            model.objective = objective_reaction
            rxn_grate[objective_reaction] = 0
        else:
            add_rxn(model,objective_reaction,target_product + ' --> ') #加入目标函数反应
            model.objective = objective_reaction
            rxn_grate[objective_reaction] = 0

        rxn_grate[objective_reaction] = 0
        solution = model.optimize()
        print(target_product,solution.fluxes[objective_reaction])

        constraint_value = 1e-6
        if solution.fluxes[objective_reaction] > 1e-6:
            met_wuxian.append(met.id)
            grate_delete_rxn = find_infinite_rxn(model,objective_reaction,constraint_value,rxn_grate)
            element_met = locating_rxn(model,check_model,objective_reaction,constraint_value,grate_delete_rxn,outputfile,unbalance_rxn)
            element_met +=element_met
    return element_met

def mac_rxn(model,rid):
    rxn = model.reactions.get_by_id(rid)
    for m in rxn.metabolites:
        if 'R' in str(m.formula) or 'X' in str(m.formula): 
            return True
        met_element_dic = m.elements
        if 'C' in met_element_dic.keys():
            if met_element_dic['C'] >= 70:
                return True
    mac_met = ['peptide','peptid','protein','PROTEIN','ppeptido','rna','DNA']
    met = [m.id.split('_')[0] for m in rxn.metabolites]
    if len(set(met) & set(mac_met)) >= 1:
        return True
    
# biomass reactions
def biomass_rxn(model,rid):
    rxn = model.reactions.get_by_id(rid)
    if re.search(r'Biomass', rxn.name) or re.search(r'biomass', rxn.name) or re.search(r'BIOMASS', rxn.id) or re.search(r'biomass', rxn.id):
        return True

def modify_delete_rxn(model,rxn_list,balance_rxn):
    cor_unri = []
    cor_unrr = []
    cor_balri = []
    cor_balrr = []
    cor_type = []

    del_ri = []
    del_rr = []
    del_type = []
    mac_bioamss = []

    for ri in rxn_list:
        rxn = model.reactions.get_by_id(ri)
        check = rxn.check_mass_balance()
        if mac_rxn(model,ri):
            del_ri.append(ri)
            del_rr.append(rxn.reaction)
            del_type.append('macromolecules')
            mac_bioamss.append(ri)
            model.remove_reactions([rxn])
            
        elif biomass_rxn(model,ri):
            del_ri.append(ri)
            del_rr.append(rxn.reaction)
            del_type.append('biomass')
            mac_bioamss.append(ri)
            model.remove_reactions([rxn])
    correct_rxn = []
    for ri in rxn_list:
        if ri not in mac_bioamss:
            rxn = model.reactions.get_by_id(ri)
            bal_cof = unbalance_coefficient(model,ri,balance_rxn)
            bal_cofactor = unbalance_miss_cofactor(model,ri,balance_rxn)
            bal_error_met = unbalance_error_met(model,ri,balance_rxn)
            
            if bal_cof:
                correct_rxn.append(ri)
                print(rxn,bal_cof,'1')
                cor_unri.append(ri)
                cor_unrr.append(rxn.reaction)
                cor_balri.append(bal_cof.id)
                cor_balrr.append(bal_cof.reaction)
                cor_type.append('coefficient')
                rxn.reaction = bal_cof.reaction
            elif bal_cofactor:
                correct_rxn.append(ri)
                print(rxn,bal_cofactor,'2')
                cor_unri.append(ri)
                cor_unrr.append(rxn.reaction)
                cor_balri.append(bal_cofactor.id)
                cor_balrr.append(bal_cofactor.reaction)
                cor_type.append('cofactor')
                rxn.reaction = bal_cofactor.reaction
            

            elif bal_error_met:
                correct_rxn.append(ri)
                cor_unri.append(ri)
                cor_unrr.append(rxn.reaction)
                cor_balri.append(bal_error_met.id)
                cor_balrr.append(bal_error_met.reaction)
                cor_type.append('metabolite error')
                rxn.reaction = bal_error_met.reaction

    for ri in rxn_list:
        if ri not in mac_bioamss and ri not in correct_rxn:
            rxn = model.reactions.get_by_id(ri)
            check = rxn.check_mass_balance()
            del_ri.append(ri)
            del_rr.append(rxn.reaction)
            del_type.append('no annotation')
            model.remove_reactions([rxn])
    return cor_unri,cor_unrr,cor_balri,cor_balrr,cor_type,del_ri,del_rr,del_type
