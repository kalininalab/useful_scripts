import ruamel.yaml
import random
import subprocess

# file locations 
PATH_to_yaml = 'config/dti/'  ## path to yaml default file
fname = "bindingDB_IC50.yaml" ## defaul yaml

# number of versions
num_versions = 60 ## how many random versions of yaml should be created

for i in range(num_versions):
    with open(PATH_to_yaml + fname, 'r') as file:
        yaml_data = ruamel.yaml.round_trip_load(file)
        # Update randomly values within the specified range
        yaml_data['model']['mlp']['dropout'] = round(random.uniform(0.01, 0.06), 2)
        # Optimizer
        yaml_data['model']['optimizer']['reduce_lr']['patience'] = random.randrange(5, 20)
        yaml_data['model']['optimizer']['reduce_lr']['factor'] = round(random.uniform(0.05, 0.3), 2)
        yaml_data['model']['optimizer']['momentum'] = round(random.uniform(0.008, 0.03), 4)
        yaml_data['model']['optimizer']['weight_decay'] = round(random.uniform(0.001, 0.01), 3)
        yaml_data['model']['optimizer']['drug_lr'] = round(random.uniform(0.0001, 0.0005), 4)
        yaml_data['model']['optimizer']['prot_lr'] = round(random.uniform(0.0001, 0.0005), 4)
        yaml_data['model']['optimizer']['lr'] = round(random.uniform(0.00035, 0.00045), 5)

    # Save the modified YAML data to a new file
    i = i+18
    new_fname = f"{fname.split('.yaml')[0]}_v{i+1}.yaml"
    with open(PATH_to_yaml + new_fname, 'w') as yaml_file:
        ruamel.yaml.round_trip_dump(yaml_data, yaml_file)

    # Run the script with the new YAML file
    command = f"python3 train.py {PATH_to_yaml}{new_fname}"
    subprocess.run(command, shell=True)
