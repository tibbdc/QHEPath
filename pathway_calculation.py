import sys
from script.bigg_function import *
if __name__ == '__main__':
    args = sys.argv[1:]

    inputpath = 'data/bigg_model_json/'
    model_respiratory_rxn_file = 'data/model_respiratory_rxn.xlsx'
    CSMN_model = 'model/CSMN.json'
    outputpath = 'results/pathway_calculation/'
    
    model_id = args[0]
    substrate = args[1]
    product = args[2]
    optimization_number = args[3] # 优化次数

    if len(args) == 5: #判断传参的个数
        outputpath = args[4]
    else:
        pass

    try:
        if not substrate.endswith('_c'): 
            if substrate[-2] == '_':
                substrate = substrate[:-1]+'c'
            if substrate[-2] != '_' and substrate[-3] == '_':
                substrate = product[:-2]+'c'
                
        if not product.endswith('_c'): 
            if product[-2] == '_':
                product = product[:-1]+'c'
            if product[-2] != '_' and product[-3] == '_':
                product = product[:-2]+'c'
        
        calculate_pathway_end(CSMN_model,model_id,substrate,product,inputpath,outputpath,model_respiratory_rxn_file,int(optimization_number))
        df_summary = get_summary_table_pathway_calculation(outputpath,substrate,product,model_id)
    except Exception as e:
        print(str(e))

