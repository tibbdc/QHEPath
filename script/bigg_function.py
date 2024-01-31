import cplex
from cplex.callbacks import MIPInfoCallback
import cobra
from cobra.flux_analysis import pfba
from cobra import Model, Reaction, Metabolite
from scipy import sparse
import numpy as np
import math
import pandas as pd
import re
import os

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

class TimeLimitCallback(MIPInfoCallback):

    def __call__(self):
        if not self.aborted and self.has_incumbent():
            gap = 100.0 * self.get_MIP_relative_gap()
            timeused = self.get_time() - self.starttime
            if timeused > self.timelimit or gap < self.acceptablegap:
                print("Good enough solution at", timeused, "sec., gap =",
                    gap, "%, quitting.")
                self.aborted = True
                self.abort()

def get_model_mat(model,model_rxn):

    """Obtaining the stoichiometric matrix of reaction coefficients in the model.

    Args:
    
    * model (cobra.Model): A Model object
    * model_rxn (list): List of reaction IDs

    return : The stoichiometric matrix of reaction coefficients in the model
    """
    mat = sparse.lil_matrix((len(model.metabolites), len(model.reactions)))
    ind = -1
    for i in model_rxn:
        ind = ind + 1
        rxn = model.reactions.get_by_id(i)
        for met, coeff in rxn.metabolites.items():
            mat[model.metabolites.index(met), ind] = coeff
    return mat

def heterogenous_onebyone(model,chassis_model,objective_reaction,N):
    """Stepwise introduction of heterologous reactions.

    Args:
    
    * model (cobra.Model): A Model object of universal model (CSMN)
    * chassis_model (cobra.Model): A Model object of host model
    * objective_reaction(str): ID of objective reaction
    * N(int): Number of optimization iterations

    """

    model.objective = objective_reaction
    v_model_rate = model.optimize().objective_value
    with chassis_model:
        chassis_model.objective = objective_reaction
        v_chassis_rate = chassis_model.optimize().objective_value
        
    # exchange and transport reactions
    exchange_rxns = []
    transport_rxns = []
    for rxn in model.reactions:
        reactants_mets = [str(m).split('_')[0] for m in rxn.reactants]
        products_mets = [str(m).split('_')[0] for m in rxn.products]

        if sorted(reactants_mets) == sorted(products_mets):
            transport_rxns.append(rxn.id)
        if len(products_mets) == 0:
            exchange_rxns.append(rxn.id) 
    n_exchange_rxns = len(exchange_rxns)
    n_transport_rxns = len(transport_rxns)
        
    rxn_all = [r.id for r in model.reactions]
    n_rxn_all = len(rxn_all)
    rxn_all_except_ex_trans = [r.id for r in model.reactions if r.id not in exchange_rxns and r.id not in transport_rxns]
    chassis_rxn1 = [r.id for r in chassis_model.reactions]
    chassis_rxn = [r.id for r in model.reactions if r.id in chassis_rxn1 and r.id not in exchange_rxns and r.id not in transport_rxns]
    n_chassis_rxn = len(chassis_rxn)
    heterogenous_rxn = [ri for ri in rxn_all_except_ex_trans if ri not in chassis_rxn]
    n_heterogenous_rxn = len(heterogenous_rxn)
    rxn_mat = exchange_rxns + transport_rxns + chassis_rxn + heterogenous_rxn
    
    # construction of A

    # S*v = 0
    mat_s = get_model_mat(model,rxn_mat)
    n_met = len(model.metabolites) 
    A = sparse.hstack([mat_s, sparse.lil_matrix((n_met, n_chassis_rxn + n_heterogenous_rxn))])
    
    #Y*lb-v <=0
    v_lb = sparse.lil_matrix(np.diag([-1] * n_rxn_all)) 
    for i in range(n_exchange_rxns + n_transport_rxns):
        v_lb[(i,i)] = 1

    y_lbr1 = sparse.lil_matrix(np.full((n_exchange_rxns+n_transport_rxns,n_chassis_rxn + n_heterogenous_rxn),0))
    y_lbr2 = sparse.lil_matrix(np.diag([-1000] * (n_chassis_rxn + n_heterogenous_rxn)))
    y_lbr = sparse.vstack([y_lbr1, y_lbr2])
    ylb = sparse.hstack([v_lb, y_lbr])
    A = sparse.vstack([A, ylb])

    # -Y*lb+v <=0
    v_ub = sparse.lil_matrix(np.diag([1] * n_rxn_all))
    y_ubr1 = sparse.lil_matrix(np.full((n_exchange_rxns+n_transport_rxns,n_chassis_rxn + n_heterogenous_rxn),0))
    y_ubr2 = sparse.lil_matrix(np.diag([-1000] * (n_chassis_rxn + n_heterogenous_rxn)))
    y_ubr = sparse.vstack([y_ubr1, y_ubr2])
    yub = sparse.hstack([v_ub, y_ubr])
    A = sparse.vstack([A, yub])
    
    
    sumy = sparse.lil_matrix([0]*n_rxn_all + [0]*n_chassis_rxn + [1]*n_heterogenous_rxn)
    A = sparse.vstack([A, sumy])
    
    
    lb = [model.reactions.get_by_id(ri).lower_bound if ri != objective_reaction else (v_chassis_rate+0.000001) for ri in rxn_mat]
    mylb = lb + [0]*(n_chassis_rxn + n_heterogenous_rxn) 

    ub = [model.reactions.get_by_id(ri).upper_bound if ri != objective_reaction else v_model_rate for ri in rxn_mat]
    myub = ub + [1]*(n_chassis_rxn + n_heterogenous_rxn) 

    my_vtypes = "C"*n_rxn_all + "B"*(n_chassis_rxn + n_heterogenous_rxn)

    var_name =['v_' + str(r) for r in rxn_mat] + ['y_' + str(r) for r in (chassis_rxn + heterogenous_rxn)]                                                                                  

    my_obj = []
    for r in rxn_mat:
        if r != objective_reaction:
            my_obj.append(0)
        else:
            my_obj.append(1)
    my_obj = my_obj + [0]*(n_chassis_rxn + n_heterogenous_rxn)
    
    my_rhs = [0]*n_met 
    myconsense = 'E'*n_met

    for ri1 in (exchange_rxns + transport_rxns):
        rxn = model.reactions.get_by_id(ri1)
        my_rhs.append(rxn.upper_bound)
        myconsense = myconsense +'L'
    for ri1 in (chassis_rxn + heterogenous_rxn):
        my_rhs.append(0)
        myconsense = myconsense +'L'

    for ri1 in (exchange_rxns + transport_rxns):
        rxn = model.reactions.get_by_id(ri1)
        my_rhs.append(rxn.upper_bound)
        myconsense = myconsense +'L'
    for ri1 in (chassis_rxn + heterogenous_rxn):
        my_rhs.append(0)
        myconsense = myconsense +'L'


    # Number of heterologous reactions
    my_rhs.append(N)
    myconsense = myconsense +'E'
    
    p = cplex.Cplex()
    timelim_cb = p.register_callback(TimeLimitCallback)
    timelim_cb.starttime = p.get_time()
    timelim_cb.timelimit = 1200
    timelim_cb.acceptablegap = 1
    timelim_cb.aborted = False

    p.parameters.tuning_status.time_limit
    p.objective.set_sense(p.objective.sense.maximize)

    p.variables.add(obj=my_obj, lb=mylb, ub=myub, types=my_vtypes,names=var_name)

    p.linear_constraints.add(rhs=my_rhs, senses=myconsense)
    rows = A.row.tolist()
    cols = A.col.tolist()
    vals = A.data.tolist()

    p.linear_constraints.set_coefficients(zip(rows, cols, vals))
    p.solve()
    
    heterogenous_rxn_dic = {}
    for ri in heterogenous_rxn:
        if round(p.solution.get_values('y_'+ri),1) == 1:
            heterogenous_rxn_dic[ri] = model.reactions.get_by_id(ri).reaction
            
    v_target = p.solution.get_objective_value()
    
    return heterogenous_rxn_dic,round(v_target,5)


def pfba_pathway_output(model,chassis_model,outputfile,fluxes,heterogenous_rxn):
    """Exporting pathways calculated by pFBA to a text (txt) file.

    Args:
    
    * model (cobra.Model): A Model object of universal model (CSMN)
    * chassis_model (cobra.Model): A Model object of host model
    * outputfile(str): Output file path
    * fluxes(dataframe): DataFrame returned by pFBA calculation
    * heterogenous_rxn(dict): Heterologous reactions(key:id of heterologous reactions; value:id of heterologous reaction equation)
    
    return: heterologous reactions introduced into the host model

    """

    flux = open(outputfile, 'w')
    flux.write('Reaction ID'+'\t'+'Fluxes'+'\t'+'Equation(ID)'+'\t'+'Equation(Name)'+'\t'+'Type'+'\n') # 文件表头
    rxn='none'
    heterogenous_rxn_pfba = {}
    for r,v in fluxes.iteritems():
        if r in heterogenous_rxn:
            rxn='heterologous'
            rxn_name = model.reactions.get_by_id(r).build_reaction_string(use_metabolite_names=True)
        else:
            rxn='native'
            rxn_name = chassis_model.reactions.get_by_id(r).build_reaction_string(use_metabolite_names=True)
        if abs(v)>1e-6:
            v=round(v,5)
            flux.write(r+'\t'+str(v)+'\t'+chassis_model.reactions.get_by_id(r).reaction+'\t'+rxn_name+'\t'+rxn+'\n')
            if rxn=='heterologous':
                heterogenous_rxn_pfba[r] = model.reactions.get_by_id(r).reaction
    flux.close()
    return heterogenous_rxn_pfba

def add_heterogenous_rxn_pfba(model,chassis_model,heterogenous_rxn,target_product,outputfile): 
    """Integrating heterologous reactions into the host model for pFBA calculation.

    Args:
    
    * model (cobra.Model): A Model object of universal model (CSMN)
    * chassis_model (cobra.Model): A Model object of host model
    * heterogenous_rxn(dict): Heterologous reactions(key:id of heterologous reactions; value:id of heterologous reaction equation)
    * target_product(str): ID of the target product
    * outputfile(str): Output file path
    
    return: heterologous reactions introduced into the host model and the flux of target product and the solution of pfba
    
    """
    with chassis_model:
        for ri,rr in heterogenous_rxn.items(): 
            rxn = model.reactions.get_by_id(ri)
            mc = [m.id.split('_')[-1] for m in rxn.metabolites]
            if mc.count('c') != len(mc):
                rr = rr.replace('_'+list(set(mc)-{'c'})[0],'_c')
                add_rxn(chassis_model,ri,rr)
            if mc.count('c') == len(mc):
                add_rxn(chassis_model,ri,rr)
        objective_reaction = 'DM_' + target_product
        solution = pfba(chassis_model)
        
        hetr = pfba_pathway_output(model,chassis_model,outputfile,solution.fluxes,heterogenous_rxn)
        bigg_product_rate = solution.fluxes[objective_reaction]

    return hetr,bigg_product_rate,solution

def add_heterogenous_rxn_fba(model,chassis_model,heterogenous_rxn,target_product): 
    """Integrating heterologous reactions into the host model for FBA calculation.

    Args:
    
    * model (cobra.Model): A Model object of universal model (CSMN)
    * chassis_model (cobra.Model): A Model object of host model
    * heterogenous_rxn(dict): Heterologous reactions(key:id of heterologous reactions; value:id of heterologous reaction equation)
    * target_product(str): ID of the target product
    
    return: heterologous reactions introduced into the host model and the flux of target product
    
    """
    with chassis_model:
        for ri,rr in heterogenous_rxn.items(): 
            rxn = model.reactions.get_by_id(ri)
            mc = [m.id.split('_')[-1] for m in rxn.metabolites]
            if mc.count('c') != len(mc):
                rr = rr.replace('_'+list(set(mc)-{'c'})[0],'_c')
                add_rxn(chassis_model,ri,rr)
            if mc.count('c') == len(mc):
                add_rxn(chassis_model,ri,rr)

        objective_reaction = 'DM_' + target_product    
        fluxes = chassis_model.optimize().fluxes

        bigg_product_rate = fluxes[objective_reaction]

    return bigg_product_rate,heterogenous_rxn

def exchange_rxn(model,rid):
    """Determining whether a reaction is an exchange reaction.

    Args:
    
    * model (cobra.Model): A Model object of universal model (CSMN)
    * rid (str): ID of reaction
    
    return: if the reaction is a exchange reaction return True
    
    """

    rxn = model.reactions.get_by_id(rid)
    products_mets = [str(m).split('_')[0] for m in rxn.products]
    if len(products_mets) == 0:
        return True
    
def transport_rxn(model,rid):
    """Determining whether a reaction is an transport reaction.

    Args:
    
    * model (cobra.Model): A Model object of universal model (CSMN)
    * rid (str): ID of reaction
    
    return: if the reaction is a transport reaction return True
    
    """
    rxn = model.reactions.get_by_id(rid)
    reactants_mets = [str(m).split('_')[0] for m in rxn.reactants]
    products_mets = [str(m).split('_')[0] for m in rxn.products]
    if sorted(reactants_mets) == sorted(products_mets):
        return True
    
def change_chassis_model_direction(model,chassis_model):
    """Align the reaction directions of the host model with the universal model(CSMN).

    Args:
    
    * model (cobra.Model): A Model object of universal model (CSMN)
    * chassis_model (cobra.Model): A Model object of host model
    
    """
    for re in chassis_model.reactions:
        if re.id in [r.id for r in model.reactions]:
            if not transport_rxn(model,re.id):
                rb = model.reactions.get_by_id(re.id)
                if re.bounds != rb.bounds:
                    re.bounds = rb.bounds    

def set_input_output(model, substrate, product, sk_bound, dm_bound):
    """Define reactions for substrate input and product output along with their respective bounds.

    Args:
    
    * model (cobra.Model): A Model object of universal model (CSMN)
    * substrate (str): ID of the substrate metabolite
    * product (str): ID of the product metabolite
    * sk_bound (float): Reaction bounds for substrate input reactions
    * dm_bound (float): Reaction bounds for product output reactions
    
    """
    if 'SK_' + substrate in model.reactions:
        model.reactions.get_by_id('SK_' + substrate).bounds=(-sk_bound,-sk_bound)
    else:
        if substrate not in model.metabolites:
            add_rxn(model,'SK_'+ substrate,substrate + ' <=> ')
            model.reactions.get_by_id('SK_' + substrate).bounds=(-sk_bound,-sk_bound)
        else:
            subm = model.metabolites.get_by_id(substrate)
            model.add_boundary(subm, type='sink',lb=-sk_bound,ub=-sk_bound)
            
    if 'DM_' + product in model.reactions:
        rxn = model.reactions.get_by_id('DM_' + product)
        rxn.reaction = product + ' --> '
        rxn.bounds = (0,dm_bound)
        model.objective = rxn
    else:
        if product not in model.metabolites:
            add_rxn(model,'DM_'+ product,product + ' --> ')
            model.reactions.get_by_id('DM_'+ product).bounds=(0,dm_bound)
            model.objective = 'DM_'+ product 
        else:
            prom = model.metabolites.get_by_id(product)
            demand = model.add_boundary(prom, type='demand', ub=dm_bound)
            model.objective = demand  

# 设置模型中除co2外，其它的含碳代谢物交换反应的边界为0
def metabolite_carbon(model):
    """Set carbon source inputs, excluding CO2, to zero.

    Args:
    
    * model (cobra.Model): A Model object of universal model (CSMN)
    
    """
    for rxn in model.reactions:
        if exchange_rxn(model,rxn.id):
            for met in rxn.metabolites:

                # Determine if it is a carbon-containing metabolite
                if re.search("C[A-Z]",str(met.formula)) or re.search("C[\d]",str(met.formula)): 
                    if met.id.split('_')[0] != 'co2' and rxn.lower_bound < 0:
                        rxn.bounds = (0,1000)

def off_optimiaze_rxn(model,rate_df):
    synthesis_rxn = []
    optimize_rxn = []
    optimize_rxn_all = []
    synthesis_rxn_all = []
    if list(rate_df['chassis_rate'])[0] != 0 and list(rate_df['het_rxn'])[0].keys() != 0: # 本源产品
        for rdic in list(rate_df['het_rxn']):
            for ri in rdic.keys():
                optimize_rxn.append(ri)
    
    if list(rate_df['chassis_rate'])[0] == 0:  # 异源产品
        for ri in list(rate_df['het_rxn'])[0]:
            synthesis_rxn.append(ri)
        
        for ri in synthesis_rxn: # 异源产品 找出synthesis_rxn包含所有分室的反应
            rxn = model.reactions.get_by_id(ri)
            hr_met1 = sorted(list(set([str(m).split('_')[0] for m in rxn.metabolites]) - {'h','h2o'}))
            for rxn2 in model.reactions:
                hr_met2 = sorted(list(set([str(m).split('_')[0] for m in rxn2.metabolites]) - {'h','h2o'}))
                if hr_met1 == hr_met2:
                    synthesis_rxn_all.append(rxn2.id)

        for rdic in list(rate_df['het_rxn'])[1:]:
            for ri in rdic.keys():
                if ri not in synthesis_rxn_all:
                    optimize_rxn.append(ri)
    for ri in optimize_rxn:   
        rxn = model.reactions.get_by_id(ri)
        hr_met1 = sorted(list(set([str(m).split('_')[0] for m in rxn.metabolites]) - {'h','h2o'}))
        for rxn2 in model.reactions:
            hr_met2 = sorted(list(set([str(m).split('_')[0] for m in rxn2.metabolites]) - {'h','h2o'}))
            if hr_met1 == hr_met2:
                optimize_rxn_all.append(rxn2.id)
    return synthesis_rxn_all,list(set(optimize_rxn_all))

def merge_respiratory(inputpath,model_id,model_respiratory_rxn_file): # inputpath 底盘模型存储路径
    """Merge the respiratory chain reactions in the host model.

    Args:
    
    * inputpath (str): The file path of the host model
    * model_id (str): The name of the host model
    * model_respiratory_rxn_file (xlsx): Excel file containing the merging of respiratory chain reactions for various types.
    
    return: A Model object of host model

    """

    chassis_model = cobra.io.load_json_model(inputpath+ model_id+'.json')

    # Proton-driven reaction in the universal model
    e_atp = ['ATPS4r','ATPS4rpp','ATPS10rpp','ATPS4m','ATPS4m_cho','ATPS3m','ATPSum','ATPS4rpp_1','ATPSu','ATPS_h','ATPS10_3','ATPS3v','ATPS3g','ATPS3r','ATPSh','ATPSu_1','ATPS4r_1','ATPS4mi']

    for rxn in chassis_model.reactions:
        if rxn.id in e_atp:
            rxn.bounds = (0,0)

    respiratory_df = pd.read_excel(model_respiratory_rxn_file,index_col='model_type')
    for mid in respiratory_df.index:
        if model_id in respiratory_df.loc[mid,'model_name']:

            for ri,rr in eval(respiratory_df.loc[mid,'respiratory_rxn']).items():
                add_rxn(chassis_model,ri,rr)

    add_rxn(chassis_model,'ADD_h','h_c <=> ') 
    return chassis_model

def product_pathway_onebyone(CSMN_model,substrate,target_product,model_id,inputpath,model_respiratory_rxn_file,outputpath):
    """Stepwise introduction of heterologous reactions for substrates and products in a host.

    Args:
    
    * CSMN_model (json): The json file of the universal model (CSMN)
    * substrate (str): ID of the substrate metabolite
    * product (str): ID of the product metabolite
    * model_id (str): The name of the host model
    * inputpath (str): The file path of the host model
    * model_respiratory_rxn_file (xlsx): Excel file containing the merging of respiratory chain reactions for various types.
    * outputpath (str): The file path of the product pathway

    return: list of threonine reaction list; list of synthesis reaction; list of synthesis reaction; list of optimization reaction
            dataframe of result; A Model object of universal model (CSMN); the number of pathways
    
    """

    model = cobra.io.load_json_model(CSMN_model)

    # 底盘模型合并呼吸链
    chassis_model = merge_respiratory(inputpath,model_id,model_respiratory_rxn_file)

    # 设置模型中除co2外，其它的含碳代谢物交换反应的边界为0
    metabolite_carbon(model)
    metabolite_carbon(chassis_model)

    change_chassis_model_direction(model,chassis_model) # 底盘模型和复合模型方向一致
    
    thr_bypass_rxn = ['THRA','THRD','R01465m','THRA_1','THRA_m','PPM2','PAI2I','PPM2_2'] # 苏氨酸循环中的反应
    with model:
        product_l = []
        chassis_rate_l = []
        asgem_rate_l = []
        het_num_l = []
        het_rxn_l = []
        thr_rxn_model = []

        set_input_output(model,substrate,target_product, 10, 1000) # 设置sink demand反应并设置通量
        set_input_output(chassis_model,substrate,target_product, 10, 1000) # 设置sink demand反应并设置通量
        objective_reaction = 'DM_' + target_product
        solution = model.optimize()
        v_model_rate = solution.objective_value
        print(target_product,v_model_rate)

        # 计算底盘模型产品得率
        chassis_solution = chassis_model.optimize()
        chassis_rate = chassis_solution.objective_value
        
        print(target_product,'chassis',chassis_rate)
        
        nr_dic = {}
        syn_rate,syn_rxn,syn_rxn_host,opt_rate,opt_rxn,opt_rxn_host = calculate_syn_opt(model,chassis_model,target_product,objective_reaction,chassis_rate,v_model_rate)
        syn_num = len(syn_rxn.keys())
        opt_num = len(opt_rxn.keys())

        nr_dic[syn_num] = syn_rxn
        nr_dic[opt_num] = opt_rxn
        print('syn_opt',syn_num,opt_num)
        
        path_num = 0
        if chassis_rate != 0 and opt_num!= 0:
            het_num = 0
            het_rxn = {}
            outputfile = outputpath+target_product+'-chassis.tsv'
            hetr_pfba,bigg_rate,pfba_solution = add_heterogenous_rxn_pfba(model,chassis_model,het_rxn,target_product,outputfile)
            product_l.append(target_product)
            chassis_rate_l.append(round(chassis_rate,5))
            asgem_rate_l.append(round(bigg_rate,5))
            het_num_l.append(het_num)
            het_rxn_l.append(het_rxn)
        
        if opt_num==0:
            het_num = 0
            het_rxn = {}
            outputfile = outputpath+target_product+'-chassis.tsv'
            hetr_pfba,bigg_rate,pfba_solution = add_heterogenous_rxn_pfba(model,chassis_model,opt_rxn,target_product,outputfile)
            product_l.append(target_product)
            chassis_rate_l.append(round(chassis_rate,5))
            asgem_rate_l.append(round(bigg_rate,5))
            het_num_l.append(het_num)
            het_rxn_l.append(het_rxn)

            # 找出是否有苏氨酸循环的反应
            for rid,v in pfba_solution.fluxes.iteritems():
                if rid in thr_bypass_rxn:
                    thr_rxn_model.append(rid)

        elif opt_num!=0:
            
            for n in range(syn_num+1,opt_num):
                try:
                    heterogenous_rxn_dic,v_product = heterogenous_onebyone(model,chassis_model,objective_reaction,n)
                    if v_product != 0:
                        nr_dic[n] = heterogenous_rxn_dic
                        print(target_product,n)
                except:
                    continue
            nv_dic = {}
            for i, heterogenous_rxn_dic in nr_dic.items():
                rate,heterogenous_rxn = add_heterogenous_rxn_fba(model,chassis_model,heterogenous_rxn_dic,target_product)
                if math.isnan(rate):
                    nv_dic[i] = 0
                else:
                    nv_dic[i] = round(rate,5)
            print(target_product,nv_dic)

            n_list = []
            if len(list(nv_dic.keys())) != 0:
                n_list.append(min(list(nv_dic.keys())))
                for n in nv_dic.keys():
                    if n <= opt_num-1:
                        if nv_dic[n+1] > nv_dic[n]:
                            n_list.append(n+1)

            # 将得率有提高的反应逐步加入到底盘模型输出途径
            path_num = len(n_list)
            m = 1
            for i in n_list:
                heterogenous_rxn_dic = nr_dic[i]
                outputfile = outputpath+target_product+'-'+str(m)+'.tsv'
                hetr_pfba,bigg_rate,pfba_solution = add_heterogenous_rxn_pfba(model,chassis_model,heterogenous_rxn_dic,target_product,outputfile)
                m = m + 1
                # 找出是否有苏氨酸循环的反应
                for rid,v in pfba_solution.fluxes.iteritems():
                    if rid in thr_bypass_rxn:
                        thr_rxn_model.append(rid)
                if len(list(hetr_pfba.keys())) !=  0:
                    het_num = len(list(hetr_pfba.keys()))
                    het_rxn = hetr_pfba
                    product_l.append(target_product)
                    chassis_rate_l.append(round(chassis_rate,5))
                    asgem_rate_l.append(round(bigg_rate,5))
                    het_num_l.append(het_num)
                    het_rxn_l.append(het_rxn)
        #将结果生成dataframe
        rate_df = pd.DataFrame({'products':product_l,'chassis_rate':chassis_rate_l,'asgem_rate':asgem_rate_l,'het_num':het_num_l,'het_rxn':het_rxn_l})
        synthesis_rxn,optimize_rxn = off_optimiaze_rxn(model,rate_df)
        return list(set(thr_rxn_model)),synthesis_rxn,optimize_rxn,rate_df,model,path_num

def pathway_iteration(i,model,substrate,target_product,model_id,inputpath,model_respiratory_rxn_file,outputpath,synthesis_rxn,optimize_rxn):
    
    chassis_model = merge_respiratory(inputpath,model_id,model_respiratory_rxn_file)

    # 关闭苏氨酸和优化反应
    thr_bypass_rxn = ['THRA','THRD','R01465m','THRA_1','THRA_m','PPM2','PAI2I','PPM2_2'] # 苏氨酸循环中的反应
    # 关闭苏氨酸循环中的反应

    for ri in thr_bypass_rxn:
        model.reactions.get_by_id(ri).bounds = (0,0)
    for ri in optimize_rxn:
            model.reactions.get_by_id(ri).bounds = (0,0)

    # 设置模型中除co2外，其它的含碳代谢物交换反应的边界为0
    metabolite_carbon(model)
    metabolite_carbon(chassis_model)
    change_chassis_model_direction(model,chassis_model) # 底盘模型和复合模型方向一致

    with model:
        set_input_output(model,substrate,target_product, 10, 1000) # 设置sink demand反应并设置通量
        set_input_output(chassis_model,substrate,target_product, 10, 1000) # 设置sink demand反应并设置通量
        objective_reaction = 'DM_' + target_product
        solution = model.optimize()
        v_model_rate = solution.objective_value

        print(target_product,v_model_rate)

        # 计算底盘模型产品得率
        chassis_solution = chassis_model.optimize()
        chassis_rate = chassis_solution.objective_value
        print(target_product,'chassis',chassis_rate)
        
        if v_model_rate == 0:
            print('no solution')
            optimize_rxn_all = []
            return list(set(optimize_rxn_all)),model
        else:
            # 得率最优时，最小异源反应步数
            opt_rxn = min_heterogenous_opt(model,chassis_model,objective_reaction,v_model_rate) # 本源优化

            # 将异源反应加入到底盘中pfba计算，并输出途径
            if len(list(opt_rxn.keys())) != 0 and sorted(list(opt_rxn.keys())) != sorted(synthesis_rxn):
                outputfile = outputpath+target_product+'-'+str(i)+'.tsv'
                hetr_pfba,bigg_rate,pfba_solution = add_heterogenous_rxn_pfba(model,chassis_model,opt_rxn,target_product,outputfile)
                optimize_rxn = [ri for ri in list(hetr_pfba.keys()) if ri not in synthesis_rxn]
                optimize_rxn_all = []
                for ri in optimize_rxn:   
                    rxn = model.reactions.get_by_id(ri)
                    hr_met1 = sorted(list(set([str(m).split('_')[0] for m in rxn.metabolites]) - {'h','h2o'}))
                    for rxn2 in model.reactions:
                        hr_met2 = sorted(list(set([str(m).split('_')[0] for m in rxn2.metabolites]) - {'h','h2o'}))
                        if hr_met1 == hr_met2:
                            optimize_rxn_all.append(rxn2.id)
                return  list(set(optimize_rxn_all)),model
            else:
                optimize_rxn_all = []
                return  list(set(optimize_rxn_all)),model
            
# 判断复合模型是否有解
def bigg_product_rate(CSMN_model,substrate,target_product):

    model = cobra.io.load_json_model(CSMN_model)

    set_input_output(model,substrate,target_product, 1, 1000) # 设置sink demand反应并设置通量
    solution = model.optimize()
    v_model_rate = solution.objective_value
    if v_model_rate == 0:
        return False
    else:
        return True

def mkdir(outputpath,substrate,product,model_id):
    outfolder = outputpath + substrate + '-'+product +'-'+model_id+'/'
    # Create a folder to save the pathway files
    path = outfolder.strip()  # Remove the space
    path = outfolder.rstrip("\\")  # Remove the '\' sign
    isExists = os.path.exists(outfolder)  # Determine whether the path exists
    if not isExists:
        os.makedirs(outfolder)
    #     # print(path + 'Created successfully')
        return outfolder
    else:
    #     # print(path + 'Directory already exists')
        return outfolder

def get_summary_table_pathway_calculation(outputpath,substrate,product,model_id):
    outfolder = outputpath + substrate + '-' + product +'-'+model_id+'/'
    path_id = []
    het_num = []
    product_yield = []
    pathway_file_name = []
    count = 1
    for file_name in os.listdir(outputpath + substrate + '-'+ product + '-'+model_id+'/'):
        pathway_df = pd.read_table(outfolder+file_name)
        pathway_df.columns = ['id','fluxes','equ','equ_name','type']
        pathway_df.set_index('id',inplace=True)
        
        if 'DM_'+product in list(pathway_df.index):
            hetrn = [k for k in pathway_df.index if pathway_df.at[k,'type'] == 'heterologous']
            path_id.append(count)
            het_num.append(len(hetrn))
            product_yield.append(round(pathway_df.at['DM_'+product,'fluxes'],4))
            pathway_file_name.append(file_name)
            count += 1
    
    model_name = [model_id] * len(path_id)
    met_df = pd.read_table('data/metabolite_id_name.tsv',index_col='metabolite_id')
    # if str(substrate) in met_df.index.to_list:
    sub_name = met_df.loc[substrate,'metabolite_name']
    # if str(product) in met_df.index.to_list:
    pro_name = met_df.loc[product,'metabolite_name']

    substrate_list = [sub_name] *len(path_id)
    product_list = [pro_name] *len(path_id)
    df_summary = pd.DataFrame({'Path ID':path_id,'Model':model_name,'Substrate':substrate_list,'Product':product_list,'No.of heterologous reactions':het_num,'Yield(mol/mol)':product_yield,'file_name':pathway_file_name})
    
    # Group the DataFrame by 'file_name' and count occurrences of zero 'No.of heterologous reactions'
    counts = df_summary.groupby('file_name')['No.of heterologous reactions'].apply(lambda x: (x == 0).sum())

    # Get a list of 'file_name' that have two or more zero 'No.of heterologous reactions'
    remove_files = counts[counts >= 2].index.tolist()

    # Filter and remove the rows that meet the conditions
    df_summary = df_summary[~((df_summary['No.of heterologous reactions'] == 0) & ~df_summary['file_name'].str.contains('chassis'))]
    df_summary = df_summary[~df_summary['file_name'].isin(remove_files)]
    
    outfolder = outputpath + substrate + '-'+ product +'-'+model_id+'/'
    # Create a folder to save the pathway files
    path = outfolder.strip()  # Remove the space
    path = outfolder.rstrip("\\")  # Remove the '\' sign
    isExists = os.path.exists(outfolder)  # Determine whether the path exists
    if not isExists:
        os.makedirs(outfolder)
    else:
        pass
    df_summary.to_csv(outfolder+'pathway_calculation_summary.tsv',sep='\t',index=False)
    return df_summary

def min_heterogenous_syn(model,chassis_model,objective_reaction,v_optimize):
    """The minimum set of heterologous reactions required to enable product synthesis.

    Args:
    
    * model (cobra.Model): A Model object of universal model (CSMN)
    * chassis_model (cobra.Model): A Model object of host model
    * objective_reaction(str): ID of objective reaction
    * v_optimize(float): The flux calculated by the universal model(CSMN)

    return: The minimum set of heterologous reactions

    """
    
    exchange_rxns = []
    transport_rxns = []
    for rxn in model.reactions:
        reactants_mets = [str(m).split('_')[0] for m in rxn.reactants]
        products_mets = [str(m).split('_')[0] for m in rxn.products]
        
        if sorted(reactants_mets) == sorted(products_mets):
            transport_rxns.append(rxn.id)
        if len(products_mets) == 0:
            exchange_rxns.append(rxn.id) 
    
    # 异源反应（除底盘之外的反应、无转运交换）chassis_model
    rxn_a1 = [r.id for r in model.reactions]
    chassis_rxn1 = [r.id for r in chassis_model.reactions]
    chassis_rxn = [r.id for r in model.reactions if r.id in chassis_rxn1]
    heterogenous_rxn_important = list(set(rxn_a1) - set(chassis_rxn) -set(exchange_rxns) -set(transport_rxns))
    heterogenous_ex_trans = list(set(rxn_a1) -set(chassis_rxn) -set(heterogenous_rxn_important))
    rxn_all = chassis_rxn + heterogenous_ex_trans + heterogenous_rxn_important

    # construction of A
    # S*v = 0
    model_rxn = rxn_all
    mat_s = get_model_mat(model,model_rxn)
    n_met = len(model.metabolites) 
    n_heterogenous_rxn_important = len(heterogenous_rxn_important)
    n_heterogenous_ex_trans = len(heterogenous_ex_trans)
    n_chassis_rxn = len(chassis_rxn)
    n_rxn_all = len(rxn_all)

    A = sparse.hstack([mat_s, sparse.lil_matrix((n_met, n_heterogenous_rxn_important))])

    #Y*lb-v <=0
    v_lb = sparse.lil_matrix(np.diag([-1] * n_rxn_all)) # np.diag([-1]*4) s生成对角矩阵
    for i in range(n_chassis_rxn + n_heterogenous_ex_trans):
        v_lb[(i,i)] = 1
        
    y_lbr1 = sparse.lil_matrix(np.full((n_rxn_all-n_heterogenous_rxn_important,n_heterogenous_rxn_important),0))
    y_lbr2 = sparse.lil_matrix(np.diag([-1000] * n_heterogenous_rxn_important))
    y_lbr = sparse.vstack([y_lbr1, y_lbr2])
    ylb = sparse.hstack([v_lb, y_lbr])
    A = sparse.vstack([A, ylb])

    # -Y*lb+v <=0
    v_ub = sparse.lil_matrix(np.diag([1] * n_rxn_all))
    y_ubr1 = sparse.lil_matrix(np.full((n_rxn_all-n_heterogenous_rxn_important,n_heterogenous_rxn_important),0))
    y_ubr2 = sparse.lil_matrix(np.diag([-1000] * n_heterogenous_rxn_important))
    y_ubr = sparse.vstack([y_ubr1, y_ubr2])
    yub = sparse.hstack([v_ub, y_ubr])
    A = sparse.vstack([A, yub])

    # construction of 

    lb = [model.reactions.get_by_id(ri).lower_bound if ri != objective_reaction else (0.01*(v_optimize)) for ri in rxn_all]
    mylb = lb + [0]*n_heterogenous_rxn_important 

    ub = [model.reactions.get_by_id(ri).upper_bound if ri != objective_reaction else (v_optimize-0.000001) for ri in rxn_all]
    myub = ub + [1]*n_heterogenous_rxn_important

    my_vtypes = "C"*n_rxn_all + "B"*n_heterogenous_rxn_important

    var_name =['v_' + str(r) for r in rxn_all] + ['y_' + str(r) for r in heterogenous_rxn_important]                                                                                  
    my_obj = [0]*n_rxn_all + [1]*n_heterogenous_rxn_important

    my_rhs = [0]*n_met # =
    myconsense = 'E'*n_met

    for ri1 in chassis_rxn + heterogenous_ex_trans:
        rxn = model.reactions.get_by_id(ri1)
        my_rhs.append(rxn.upper_bound)
        myconsense = myconsense +'L'
    for ri1 in heterogenous_rxn_important:
        my_rhs.append(0)
        myconsense = myconsense +'L'
        
    for ri1 in chassis_rxn + heterogenous_ex_trans:
        rxn = model.reactions.get_by_id(ri1)
        my_rhs.append(rxn.upper_bound)
        myconsense = myconsense +'L'
    for ri1 in heterogenous_rxn_important:
        my_rhs.append(0)
        myconsense = myconsense +'L'
    p = cplex.Cplex()
    p.objective.set_sense(p.objective.sense.minimize)

    p.variables.add(obj=my_obj, lb=mylb, ub=myub, types=my_vtypes,names=var_name)

    p.linear_constraints.add(rhs=my_rhs, senses=myconsense)
    rows = A.row.tolist()
    cols = A.col.tolist()
    vals = A.data.tolist()

    p.linear_constraints.set_coefficients(zip(rows, cols, vals))
    p.solve()
    
    add_heterogenous_rxn = {}
    for ri in heterogenous_rxn_important:
        if round(p.solution.get_values('y_'+ri),1) == 1:
            add_heterogenous_rxn[ri] = model.reactions.get_by_id(ri).reaction
    return add_heterogenous_rxn

def min_heterogenous_opt(model,chassis_model,objective_reaction,v_optimize):
    """The minimum set of heterologous reactions required for product optimization.

    Args:
    
    * model (cobra.Model): A Model object of universal model (CSMN)
    * chassis_model (cobra.Model): A Model object of host model
    * objective_reaction(str): ID of objective reaction
    * v_optimize(float): The flux calculated by the universal model(CSMN)

    return: The minimum set of heterologous reactions

    """

    exchange_rxns = []
    transport_rxns = []
    for rxn in model.reactions:
        reactants_mets = [str(m).split('_')[0] for m in rxn.reactants]
        products_mets = [str(m).split('_')[0] for m in rxn.products]
        
        if sorted(reactants_mets) == sorted(products_mets):
            transport_rxns.append(rxn.id)
        if len(products_mets) == 0:
            exchange_rxns.append(rxn.id) 
    
    # 异源反应（除底盘之外的反应、无转运交换）chassis_model
    rxn_a1 = [r.id for r in model.reactions]
    chassis_rxn1 = [r.id for r in chassis_model.reactions]
    chassis_rxn = [r.id for r in model.reactions if r.id in chassis_rxn1]

    heterogenous_rxn_important = list(set(rxn_a1) - set(chassis_rxn))

    rxn_all = chassis_rxn + heterogenous_rxn_important

    # construction of A
    # S*v = 0
    model_rxn = rxn_all
    mat_s = get_model_mat(model,model_rxn)
    n_met = len(model.metabolites) 
    n_heterogenous_rxn_important = len(heterogenous_rxn_important)
    n_chassis_rxn = len(chassis_rxn)
    n_rxn_all = len(rxn_all)

    A = sparse.hstack([mat_s, sparse.lil_matrix((n_met, n_heterogenous_rxn_important))])

    #Y*lb-v <=0
    v_lb = sparse.lil_matrix(np.diag([-1] * n_rxn_all)) # np.diag([-1]*4) s生成对角矩阵
    for i in range(n_chassis_rxn):
        v_lb[(i,i)] = 1
        
    y_lbr1 = sparse.lil_matrix(np.full((n_rxn_all-n_heterogenous_rxn_important,n_heterogenous_rxn_important),0))
    y_lbr2 = sparse.lil_matrix(np.diag([-1000] * n_heterogenous_rxn_important))
    y_lbr = sparse.vstack([y_lbr1, y_lbr2])
    ylb = sparse.hstack([v_lb, y_lbr])
    A = sparse.vstack([A, ylb])

    # -Y*lb+v <=0
    v_ub = sparse.lil_matrix(np.diag([1] * n_rxn_all))
    y_ubr1 = sparse.lil_matrix(np.full((n_rxn_all-n_heterogenous_rxn_important,n_heterogenous_rxn_important),0))
    y_ubr2 = sparse.lil_matrix(np.diag([-1000] * n_heterogenous_rxn_important))
    y_ubr = sparse.vstack([y_ubr1, y_ubr2])
    yub = sparse.hstack([v_ub, y_ubr])
    A = sparse.vstack([A, yub])

    # construction of 
    lb = [model.reactions.get_by_id(ri).lower_bound if ri != objective_reaction else v_optimize for ri in rxn_all]
    mylb = lb + [0]*n_heterogenous_rxn_important 

    ub = [model.reactions.get_by_id(ri).upper_bound if ri != objective_reaction else v_optimize for ri in rxn_all]
    myub = ub + [1]*n_heterogenous_rxn_important

    my_vtypes = "C"*n_rxn_all + "B"*n_heterogenous_rxn_important

    var_name =['v_' + str(r) for r in rxn_all] + ['y_' + str(r) for r in heterogenous_rxn_important]                                                                                  
    my_obj = [0]*n_rxn_all + [1]*n_heterogenous_rxn_important

    my_rhs = [0]*n_met # =
    myconsense = 'E'*n_met

    for ri1 in chassis_rxn:
        rxn = model.reactions.get_by_id(ri1)
        my_rhs.append(rxn.upper_bound)
        myconsense = myconsense +'L'
    for ri1 in heterogenous_rxn_important:
        my_rhs.append(0)
        myconsense = myconsense +'L'
        
    for ri1 in chassis_rxn:
        rxn = model.reactions.get_by_id(ri1)
        my_rhs.append(rxn.upper_bound)
        myconsense = myconsense +'L'
    for ri1 in heterogenous_rxn_important:
        my_rhs.append(0)
        myconsense = myconsense +'L'
    p = cplex.Cplex()

    timelim_cb = p.register_callback(TimeLimitCallback)
    timelim_cb.starttime = p.get_time()
    timelim_cb.timelimit = 600
    timelim_cb.acceptablegap = 1
    timelim_cb.aborted = False
    p.parameters.tuning_status.time_limit

    p.objective.set_sense(p.objective.sense.minimize)

    p.variables.add(obj=my_obj, lb=mylb, ub=myub, types=my_vtypes,names=var_name)

    p.linear_constraints.add(rhs=my_rhs, senses=myconsense)
    rows = A.row.tolist()
    cols = A.col.tolist()
    vals = A.data.tolist()

    p.linear_constraints.set_coefficients(zip(rows, cols, vals))
    p.solve()
    
    add_heterogenous_rxn = {}
    for ri in heterogenous_rxn_important:
        if round(p.solution.get_values('y_'+ri),1) == 1:
            add_heterogenous_rxn[ri] = model.reactions.get_by_id(ri).reaction
    return add_heterogenous_rxn

def add_heterogenous_rxn_pfba_syn_opt(model,chassis_model,heterogenous_rxn,target_product): 
    with chassis_model:
        objective_reaction = 'DM_' + target_product
        for ri,rr in heterogenous_rxn.items(): 
            add_rxn(chassis_model,ri,rr)
        chassis_model.objective=objective_reaction
        # solution = chassis_model.optimize()
        solution = pfba(chassis_model)
        add_heterogenous_rxn = {}
        bigg_product_rate = solution.fluxes[objective_reaction]
        for rid,v in solution.fluxes.iteritems():
            if abs(v) > 1e-6:
                if rid in list(heterogenous_rxn.keys()):
                    add_heterogenous_rxn[rid] = heterogenous_rxn[rid]
    return bigg_product_rate,add_heterogenous_rxn

# 计算产品的合成优化得率
def calculate_syn_opt(model,chassis_model,target_product,objective_reaction,v_chassis_rate,v_model_rate): 
    """Calculate the synthesis and optimize the yield of the product in the host model after introducing heterologous reactions.

    Args:
    
    * model (cobra.Model): A Model object of universal model (CSMN)
    * chassis_model (cobra.Model): A Model object of host model
    * target_product (str): ID of the target product
    * objective_reaction(str): ID of objective reaction
    * v_chassis_rate(float): The flux calculated by the host model
    * v_model_rate(float): The flux calculated by the universal model(CSMN)

    return: 
        syn_rate(float): flux of product synthesis
        syn_rxn(dict): The minimum set of heterologous reactions required for product synthesis
        syn_rxn_host(dict): The minimum set of heterologous reactions required for product synthesis except exchange and transport reactions
        opt_rate(float): flux of product synthesis
        opt_rxn(dict): The minimum set of heterologous reactions required for product optimization
        opt_rxn_host(dict): The minimum set of heterologous reactions required for product optimization except exchange and transport reactions
        
    """

    if v_chassis_rate > 0: 

        syn_rate = v_chassis_rate # 本源合成
        syn_rxn = {}
        syn_rxn_host = {}
        
        opt_rxn = min_heterogenous_opt(model,chassis_model,objective_reaction,v_model_rate) # 本源优化
        opt_rate, opt_rxn = add_heterogenous_rxn_pfba_syn_opt(model,chassis_model,opt_rxn,target_product)
        del_trans_rxn = [ri for ri in opt_rxn.keys() if transport_rxn(model,ri)] # 去掉交换反应
        del_exchange_rxn = [ri for ri in opt_rxn.keys() if exchange_rxn(model,ri)] # 去掉转运反应
        opt_rxn_host = dict(opt_rxn)
        for ri in del_exchange_rxn:
            opt_rxn_host.pop(ri)
        for ri in del_trans_rxn:
            opt_rxn_host.pop(ri)
        
        return round(syn_rate,5),syn_rxn,syn_rxn_host,round(opt_rate,5),opt_rxn,opt_rxn_host

    else:
        syn_rxn = min_heterogenous_syn(model,chassis_model,objective_reaction,v_model_rate) #非本源合成
        syn_rate,syn_rxn = add_heterogenous_rxn_pfba_syn_opt(model,chassis_model,syn_rxn,target_product)
        del_trans_rxn = [ri for ri in syn_rxn.keys() if transport_rxn(model,ri)] # 去掉交换反应
        del_exchange_rxn = [ri for ri in syn_rxn.keys() if exchange_rxn(model,ri)] # 去掉转运反应
        syn_rxn_host = dict(syn_rxn)
        for ri in del_exchange_rxn:
            syn_rxn_host.pop(ri)
        for ri in del_trans_rxn:
            syn_rxn_host.pop(ri)

        opt_rxn = min_heterogenous_opt(model,chassis_model,objective_reaction,v_model_rate) # 非本源优化
        opt_rate,opt_rxn = add_heterogenous_rxn_pfba_syn_opt(model,chassis_model,opt_rxn,target_product)
        del_trans_rxn = [ri for ri in opt_rxn.keys() if transport_rxn(model,ri)] # 去掉交换反应
        del_exchange_rxn = [ri for ri in opt_rxn.keys() if exchange_rxn(model,ri)] # 去掉转运反应
        opt_rxn_host = dict(opt_rxn)
        for ri in del_exchange_rxn:
            opt_rxn_host.pop(ri)
        for ri in del_trans_rxn:
            opt_rxn_host.pop(ri)
        return round(syn_rate,5),syn_rxn,syn_rxn_host,round(opt_rate,5),opt_rxn,opt_rxn_host

def calculate_pathway_end(CSMN_model,model_id,substrate,product,inputpath,outputpath,model_respiratory_rxn_file,n):
    
    outfolder = outputpath + substrate + '-' + product +'-'+model_id+'/'
    # Create a folder to save the pathway files
    path = outfolder.strip()  # Remove the space
    path = outfolder.rstrip("\\")  # Remove the '\' sign
    isExists = os.path.exists(outfolder)  # Determine whether the path exists
    if not isExists:
        os.makedirs(outfolder)
    else:
        pass

    if bigg_product_rate(CSMN_model,substrate,product): # 判断是否有解
        thr_rxn,synthesis_rxn,optimize_rxn,rate_df,model,path_num  = product_pathway_onebyone(CSMN_model,substrate,product,model_id,inputpath,model_respiratory_rxn_file,outfolder)
        print(optimize_rxn)

        # n = 6 #迭代次数
        # for i in range(path_num+1,path_num+1+n):
        for i in range(path_num+1,path_num+n):
            optimize_rxn,model = pathway_iteration(i,model,substrate,product,model_id,inputpath,model_respiratory_rxn_file,outfolder,synthesis_rxn,optimize_rxn)
            print(i,optimize_rxn)
            if len(optimize_rxn) == 0:
                break
    else:
        print('no solution')