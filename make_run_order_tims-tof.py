import pandas as pd
import re
import numpy as np
 

raw_file_dest = "D:\\Data\\Test\\You Forgot"

group_prefixes = [
    "p1.6",
    "p1.8"
    "Blank"
    ] #these goes into the raw file name
channel_prefixes = ["_CH1_","_CH2_","_CH3_"] #this too

group_methods = [
    ["channel_HCD30_1.meth", "channel_HCD30_2.meth","channel_HCD30_3.meth"], 
    ["channel_HCD30_1.meth", "channel_HCD30_2.meth","channel_HCD30_3.meth"],
    ["channel_HCD20_1.meth", "channel_HCD20_2.meth","channel_HCD20_3.meth"],
    ["channel_HCD20_1.meth", "channel_HCD20_2.meth","channel_HCD20_3.meth"],
    ["channel_HCD30_1.meth", "channel_HCD30_2.meth","channel_HCD30_3.meth"],
    ]
INJ_VOLUME = 1 #this is not really used for us

group_sample_wells = [
    ["RA3", "RA4", "RA13"],
    ["RB2", "RB5", "RB11"],
    "Hela_hcd20_wells.csv",
    "K562_hcd20_wells.csv", #if you put your wells in a csv, with a column called "Position", it can read that in
    "7@R5" #if always the same location put the number and location a string separated with an @
]

group_suffixes = [
    "",
    "",
    "",
    "",
    "numerical"    
] #numerical is a keyword that will use "_" and the rep number,
  #use "" for no suffix

OUTPUT_FILENAME = "output.csv"

Block_over_time = True
Randomize = True
Block_between_channels = True



################################### Don't change anything beyond this point ##################################
# Define functions

#Block over time
def add_blocks(table):
    num_samples = table.shape[0]
    new_table = table.copy()
    new_table["Block Width"] = 1/num_samples
    new_table["Block"] = range(0,num_samples)
    return new_table

#Randomize only
def add_randomizer(table):
    num_samples = table.shape[0]
    new_table = table.copy()
    new_table["Block Width"] = 1/num_samples
    new_table["Randomizer"] = [np.random.rand() for x in range(num_samples)]
    return new_table

def format_tables(tables, randomize, block):
    formatted_tables = []
    for each_table in tables:
        current_table = each_table.copy()
        if block:
            current_table = add_blocks(current_table)
        if randomize:
            current_table = add_randomizer(current_table)
        formatted_tables.append(current_table)
    return formatted_tables

#Reorder
def reorder_table(table):
    #Both were done
    if  "Block" in table.columns and \
        "Block Width" in table.columns and \
        "Randomizer" in table.columns:
                table["Order"] = table["Block"] * table["Block Width"] + \
                                 table["Randomizer"] * table["Block Width"]
    elif "Block Width" in table.columns and \
        "Randomizer" in table.columns:
                table["Order"] = table["Randomizer"] * table["Block Width"]
    elif "Block Width" in table.columns and \
        "Block" in table.columns:
                table["Order"] = table["Block"] * table["Block Width"]
    table = table.sort_values(by="Order").reset_index(drop = True)
    return table

#remove extra columns
def tidy_table(table):
    for each_extra_column in ["Order","Block","Block Width","Randomizer"]:
        if each_extra_column in table.columns:
            table = table.drop(each_extra_column, axis=1)
    return table

#check you did it right
if len(group_prefixes) > len(group_sample_wells):
    print("Missing sample wells")
    exit()
elif len(group_prefixes) < len(group_sample_wells):
    print("extra sample wells")
    exit()
i = 0
for each_groups_methods in group_methods:
    i = i + 1
    if len(each_groups_methods) != len(channel_prefixes):
        print("not enough methods for group {i}.")
        exit()

if len(group_prefixes) != len(group_suffixes):
    print('wrong number of suffixes\nit should match the number of prefixes, fill with "".')

tables = []
j = 0
for each_prefix in group_prefixes:
    current_wells = group_sample_wells[j]
    current_type = type(current_wells)
    num_samples = len(current_wells)
    if current_type == str and "." in current_wells: #file
        current_data_frame = pd.read_csv(current_wells)
    elif current_type == str and "@" in current_wells: #single well
        num_samples = int(current_wells.split("@")[0])
        current_wells = current_wells.split("@")[1]
        current_data_frame = pd.DataFrame(data={"Position": [current_wells for x in range(num_samples)]})
    elif current_type == list:
        current_data_frame = pd.DataFrame(data={"Position": current_wells})
    else:
        print(current_wells)
    well_container = re.sub(pattern="[0-9]+",repl="",string=current_data_frame["Position"].tolist()[0])
    if len(well_container) == 1: #buffer
        current_data_frame["Sample Type"] = "Blank"
    if len(well_container) == 2: #wellplate
        current_data_frame["Sample Type"] = "Unknown"

    current_data_frame["File Name"] = each_prefix
    current_data_frame["Path"] = raw_file_dest
    current_data_frame["Inj Vol"] = INJ_VOLUME #injection volume
    
    tables.append(current_data_frame)
    j = j + 1



if Block_between_channels:
    num_channels = len(channel_prefixes)
    #split into channels
    channels = []
    for each_channel in range(num_channels):
        channels.append([])
        for each_table in tables:
            num_rows = each_table.shape[0]
            sliced_table = each_table.iloc[[x for x in range(num_rows) \
                                            if x%num_channels == each_channel],:]
            if sliced_table.shape[0] > 0:
                channels[each_channel].append(sliced_table)

    
    
    #format and combine each column
    formatted_channels = []
    for each_channel in channels:
        formatted_tables = format_tables(each_channel, Randomize, Block_over_time)
        combine_table = pd.DataFrame()
        for each_table in formatted_tables:
            combine_table = pd.concat([combine_table,each_table])
        
        formatted_channels.append(combine_table)

    #redistribute
    current_distribution = [x.shape[0] for x in formatted_channels]
    samples_per_channel = int(sum(current_distribution)/num_channels)
    leftover_samples = sum(current_distribution)%num_channels
    
    desired_distribution = [samples_per_channel for x in range(num_channels)]
    for x in range(leftover_samples):
        desired_distribution[x] += 1

    for x in range(len(formatted_channels) -1):
        final_index = desired_distribution[x]
        formatted_channels[x+1] = pd.concat([formatted_channels[x+1],
                                             formatted_channels[x].iloc[final_index:,:]])
        formatted_channels[x] = formatted_channels[x].iloc[0:final_index,:]

    #reorder 
    sorted_channels = []
    for each_channel in formatted_channels:
        sorted_channels.append(reorder_table(each_channel))


    #intercalate samples from each column
    
    final_table = pd.DataFrame()
    for each_set in range(samples_per_channel):
        for each_channel in sorted_channels:
            current_table = each_channel.iloc[[each_set],:]
            final_table = pd.concat([final_table, current_table])
    if desired_distribution[0] > samples_per_channel:
        for each_channel in sorted_channels:
            if each_channel.shape[0] > samples_per_channel:
                current_table = each_channel.iloc[[-1],:]
                final_table = pd.concat([final_table, current_table])
    final_table = final_table.reset_index(drop=True)
    
    #rename filename and add method
    for index, row in final_table.iterrows():
        old_filename = row["File Name"]
        current_channel = index%num_channels
        current_suffix_index = group_prefixes.index(old_filename)
        current_suffix =  group_suffixes[current_suffix_index]
        #see if suffix is sample number
        if current_suffix == "numerical" and "Block" in row.keys():            
            current_suffix = "_" + str(row["Block"] + 1)  

        new_filename =  old_filename\
                + channel_prefixes[current_channel] \
                + row["Position"]\
                + current_suffix
        final_table.at[index,"File Name"] = new_filename
        final_table.at[index,"Instrument Method"] = group_methods[current_suffix_index][current_channel]

else:
    #format tables
    formatted_tables = format_tables(tables, Randomize, Block_over_time)

    #combine tables
    combine_table = pd.DataFrame()
    for each_table in formatted_tables:
        combine_table = pd.concat([combine_table,each_table])

    final_table = reorder_table(combine_table)

    #rename filename and add method
    for index, row in final_table.iterrows():
        old_filename = row["File Name"]
        current_channel = index%len(channel_prefixes)
        current_suffix_index = group_prefixes.index(old_filename)
        current_suffix =  group_suffixes[current_suffix_index]
        #see if suffix is sample number
        if current_suffix == "numerical" and "Block" in row.keys():            
            current_suffix = "_" + str(row["Block"] + 1)  

        new_filename =  old_filename\
                + channel_prefixes[current_channel] \
                + row["Position"]\
                + current_suffix
        final_table.at[index,"File Name"] = new_filename
        final_table.at[index,"Instrument Method"] = group_methods[current_suffix_index][current_channel]


#make column method and add channels
final_table = tidy_table(final_table)

#Export
# write the header
header = 'Bracket Type=4.\n'

with open(OUTPUT_FILENAME, 'w') as fp:
    fp.write(header)

# write the rest
final_table.to_csv(OUTPUT_FILENAME, header=True, mode='a', index=False)