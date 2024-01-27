import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import os
from pathlib import Path


INPUT_FILES = {"10 min":"10min_January/report.tsv","20 min":"20min_January/report.tsv","18 min":"report.tsv"}
# INPUT_FILES = {"20 min":"report.tsv"}
File_Category_Name = "Gradient Length"
FILTERS = ["10ng","16ng","1ng","200pg"]
# FILTERS = ["16ng","1ng","200pg"]
Filter_Category_Name = "Sample Amount"
WRITE_OUTPUT = True
APPFOLDER = "./"
# output_name = "ThyRTs"
output_name = "All Stdev"


def plot_RTstdev_boxplot(RTs, plot_options, username):
    
    plot_div = None
    CSV_link = None
    SVG_link = None

    if plot_options["log10"]:
        RTs= RTs[RTs["Stdev"]>0]

    if not plot_options["ylimits"] or plot_options["ylimits"] == "[]" or \
            not isinstance(plot_options["ylimits"], list):
        ylimits = None
    else:
        ylimits = plot_options["ylimits"]

    #Remove outliers
    if plot_options["outliers"] == "remove":
        q_low = RTs["Stdev"].quantile(0.01)
        q_hi  = RTs["Stdev"].quantile(0.99)

        RTs = RTs[(RTs["Stdev"] < q_hi) & (RTs["Stdev"] > q_low)]
        plot_options["outliers"] = False
        
    RT_summary = RTs.groupby(["Group Name"]).agg(
        meds = ("Stdev","median")).reset_index()
    # median label
    if plot_options["median label"] == "True" or \
            plot_options["median label"] == True:
        total_labels = [{"x": x, "y": total*1.15, "text": str(
            round(total,1)), "showarrow": False} for x, total in zip(
            RT_summary["Group Name"], RT_summary["meds"])]
    else:
        total_labels = []   # no median labels

    if plot_options["RT mode"] == "grouped":  
            
        #find out present categories
        categories = RTs.groupby(plot_options["Group By Color"]).first().reset_index()[plot_options["Group By Color"]].tolist()
        # create the plot
        fig_data = []
        i = 0
        for eachCategory in categories:
            fig_data.append(go.Box(name = eachCategory,
                        x=RTs.loc[RTs[plot_options["Group By Color"]]==eachCategory,plot_options["Group By X"]].tolist(),
                        y=RTs.loc[RTs[plot_options["Group By Color"]]==eachCategory,"Stdev"].tolist(),
                        fillcolor = plot_options["color"][i],
                        boxpoints=plot_options["outliers"],
                        marker=dict(opacity=0)
                        ))
            i = i + 1

        fig = go.Figure(data = fig_data)
        fig.update_layout(
            boxmode = "group",
            plot_bgcolor='white',
            paper_bgcolor='white',
            yaxis=dict(title=plot_options["Y Title"],showline=True, linewidth=1, linecolor='black'),
            xaxis=dict(title=plot_options["X Title"],showline=True, linewidth=1, linecolor='black')
            )
        
    else:
    # create the interactive plot
        fig = px.box(RTs,
                        x="Group Name",
                        y='Stdev',
                        color="Group Name",
                        color_discrete_sequence=plot_options["color"],
                        width=plot_options["width"],
                        height=plot_options["height"],
                        )

        fig.update_layout(
            yaxis=dict(title=plot_options["Y Title"],
                    range=ylimits, showline=True, linewidth=1, linecolor='black'),
            font=plot_options["font"],
            xaxis=dict(title=plot_options["X Title"],showline=True, linewidth=1, linecolor='black'),
            showlegend=True,
            annotations=total_labels,
            boxmode = "group",
            plot_bgcolor='white',
            paper_bgcolor='white',
        )
    if WRITE_OUTPUT:        
        # create the file for donwnload
        img_dir = os.path.join(APPFOLDER, "images/")
        if not os.path.exists(img_dir):
            Path(img_dir).mkdir(parents=True)

        fig.write_image(os.path.join(
            img_dir, f"{username}_RT_Boxplot.svg"), format = "svg", validate = False, engine = "kaleido")
        
        # create the download CSV and its link
        data_dir = os.path.join(APPFOLDER, "csv/")
        if not os.path.exists(data_dir):
            Path(data_dir).mkdir(parents=True)
        RTs.to_csv(os.path.join(
            data_dir, f"{username}_all_RTs.csv"), index=False)
        RT_summary.to_csv(os.path.join(
            data_dir, f"{username}_RT_Summary.csv"), index=False)
        print("Downloading links...")
        


    return fig

all_names = []
i=0
for eachType in INPUT_FILES.keys():
    eachFile = INPUT_FILES[eachType]
    allData = pd.read_table(eachFile,sep="\t")
    
    for eachFilter in FILTERS:
        currentData = allData[allData["Run"].str.contains(eachFilter)]
        currentData = currentData[~currentData["Run"].str.contains("1775")]
        currentData = currentData.groupby("Precursor.Id").agg(
                                                                    # Median= ("RT","median"),
                                                                    # Mean=("RT","mean"),
                                                                    Stdev=("RT",np.std)).reset_index()
        # currentData[File_Category_Name] = eachType
        # currentData[Filter_Category_Name] = eachFilter
        # currentData["Group Name"] = eachType + " " + eachFilter
        new_names = [str(x) + " " + eachType + " " + eachFilter for x in currentData.columns.to_series()[1:]]
        
        mapping  = dict(zip(currentData.columns.to_series()[1:], new_names))
        currentData = currentData.rename(columns=mapping)
        if currentData.size >0:
            all_names.append("Stdev "+eachType + " " + eachFilter)
            if i == 0:
                allPeaks = currentData
            else:
                allPeaks = pd.merge(allPeaks,currentData,how="inner")
                # allPeaks = pd.concat([allPeaks,currentData])
        i = i + 1

# print(all_names)
# print(allPeaks)
allPeaks = allPeaks.melt(value_vars=all_names,
              id_vars="Precursor.Id",var_name="Group Name",value_name="Stdev")

# print(allPeaks)

plot_options={    
        "median label": "True",
        "X Title": "Gradient Length",
        "Y Title": "Standard Deviation in Retention Time",
        "color": ["#F39B6D","#fac8d3", "#7aa8e6","#3E6990","blue", "cyan", "black", "gray", "purple", "pink",
                  "orange", "brown", "pink", "gray", "olive", "cyan"],
        "log10": True,
        "width": 700,
        "height": 450,
        "font": dict(size=32, family="Arial black"),
        "outliers": False, #False, "all", "outliers", "suspectedoutliers", "remove"
        "ylimits": [],
        "RT mode": "ungrouped", #grouped, ungrouped
        "Group By X": "Sample Amount",
        "Group By Color": "SPE Type",
        "help for information only": \
        "median label options: True or False." \
        "color: the first few colors will be used"\
    }
fig = plot_RTstdev_boxplot(allPeaks, plot_options,username=output_name)
fig.update_traces(boxpoints=False,showlegend=False)
fig.write_image("images/"+output_name+".png", scale=10, width = 1000, height=1000) 
fig.show()
