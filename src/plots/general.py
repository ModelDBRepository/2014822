from ..helper import debug_reader as dr
import matplotlib.pyplot as plt
import numpy as np
import os
from pprint import pprint
import matplotlib.image as mpimg
from src.helper.config import load_config, load_debug_config

dpi = 300
path = "figures"

def plot_by_datestring(datestrings, filenamesuffix, protocol, allparams):
    datas = {}
    for datestring in datestrings:
        data = dr.read_protocol_bydatestring(datestring, protocol)
        params_data = [data[allparams[0]], data[allparams[1]]]
        datas[datestring.split('.')[-1]] = params_data

    for item in datas:
        print(20*'=', item)

        print(f"LTP: {datas[item][0]['epsp']:.4f}")
        print(f"LTD: {datas[item][1]['epsp']:.4f}")

    datetime = ".".join(datestrings[0].split('.')[:-1])
    config = load_debug_config(datetime)

    extension = '.png'
    # plot_hist
    hist_vmem_filename = f"{datestring}_hist_vmem{filenamesuffix}_test{extension}"
    datanames = config['datanames']['hist']
    hist_vmem_nr2b_data = {x: datas[x] for x in datas if x in datanames}
    plot_hist_vmem_nr2b(hist_vmem_nr2b_data, hist_vmem_filename, datanames, config, 'hist', savemore=f"f1.hist_vmem{extension}")

    # plot_hist
    hist_vmem_filename = f"{datestring}_hist_vmem_restored{filenamesuffix}_test{extension}"
    datanames = config['datanames']['histr']
    hist_vmem_nr2b_data = {x: datas[x] for x in datas if x in datanames}
    plot_hist_vmem_nr2b(hist_vmem_nr2b_data, hist_vmem_filename, datanames, config, 'histr', savemore=f"f2.hist_vmem_restored{extension}")

    # plot_epsp
    datanames = config['datanames']['epsp']
    epsp_experiment_filename = f"{datestring}_epsp_experiment{filenamesuffix}_test{extension}"
    epsp_experiment_data = {x: datas[x] for x in datas if x in datanames}
    plot_epsp_experiment(epsp_experiment_data, epsp_experiment_filename, datanames, config, savemore=f"f6.epsp_experiment{extension}")
    plot_epsp_experiment(epsp_experiment_data, epsp_experiment_filename, datanames, config, savemore=f"f6.ticks_epsp_experiment{extension}", epspticks=True)

    # plot_epsp nr2baicd
    configname = 'nr2baicd'
    datanames = config['datanames'][configname]
    epsp_experiment_filename = f"{datestring}_epsp_{configname}{filenamesuffix}_test{extension}"
    epsp_experiment_data = {x: datas[x] for x in datas.keys() if x in datanames}
    plot_epsp_nr2b(epsp_experiment_data, epsp_experiment_filename, datanames, config, configname, savemore=f"f4.epsp_nr2baicd{extension}", inset=True)
    plot_epsp_nr2b(epsp_experiment_data, epsp_experiment_filename, datanames, config, configname, savemore=f"f4.ticks_epsp_nr2baicd{extension}", epspticks=True, inset=True)

    # plot_epsp nr2bbeta
    configname = 'nr2bbeta'
    datanames = config['datanames'][configname]
    epsp_experiment_filename = f"{datestring}_epsp_{configname}{filenamesuffix}_test{extension}"
    epsp_experiment_data = {x: datas[x] for x in datas.keys() if x in datanames}
    plot_epsp_nr2b(epsp_experiment_data, epsp_experiment_filename, datanames, config, configname, savemore=f"f3.epsp_nr2bbeta{extension}")
    plot_epsp_nr2b(epsp_experiment_data, epsp_experiment_filename, datanames, config, configname, savemore=f"f3.ticks_epsp_nr2bbeta{extension}", epspticks=True)
    
    # plot_epsp nr2baicdbeta
    configname = 'nr2baicdbeta'
    datanames = config['datanames'][configname]
    epsp_experiment_filename = f"{datestring}_epsp_{configname}{filenamesuffix}_test{extension}"
    epsp_experiment_data = {x: datas[x] for x in datas.keys() if x in datanames}
    plot_epsp_nr2b(epsp_experiment_data, epsp_experiment_filename, datanames, config, configname, savemore=f"f5.epsp_nr2baicdbeta{extension}")
    plot_epsp_nr2b(epsp_experiment_data, epsp_experiment_filename, datanames, config, configname, savemore=f"f5.ticks_epsp_nr2baicdbeta{extension}", epspticks=True)



def plot_hist_vmem_nr2b(datas, filename, datanames, config, configlabel, savemore=None):
    labels = config['labels'][configlabel]
    colors = config['colors'][configlabel]
    colors = np.array(colors)

    title_size = 32
    tick_size = 18
    label_size = 28
    label_size2 = 32
    legend_size = 26
    linewidth = 5
    
    start_end = [
        [1300, 3000],
        [1300, 3000],
    ]

    vmem_ylims = [
        [-70, -50],
        [-66, -55],
    ]
    
    glu_ylims = [
        [-0.01e-3, 1.72e-3],
        [-0.01e-3, 0.704e-3]
    ]

    epsp_ylims = [
        [-64, -59],
        [-64, -59],
    ]

    epsp_diff = 30
    epsp_xlim = 30
    epsp_start_end = [200 - epsp_diff, -200 - epsp_diff]

    for protocol_idx, protocol_name in enumerate(['ltp', 'ltd']):
        fig, axes = plt.subplots(
            nrows=len(datanames), ncols=4,
            figsize=(52, len(datanames)*8),
            width_ratios=[5, 5, 2, 3])
        for idx, run_name in enumerate(datanames):
            data = datas[run_name]
            data = data[protocol_idx]

            plot_start = start_end[protocol_idx][0]
            plot_end = start_end[protocol_idx][1]
            
            xindexes = [i for i, x in enumerate(list(data['t'])) if x > plot_start and x < plot_end]
            epsp_xindexes = [i for i, x in enumerate(list(data['t'])) if x > epsp_start_end[0] and (x < epsp_start_end[0] + epsp_diff*2)]
            epsp_xindexes2 = [i for i, x in enumerate(list(data['t'])) if (x > (epsp_start_end[1] + data['t'][-1])) and (x < (epsp_start_end[1] + epsp_diff*2 + data['t'][-1]))]

            epsp_xindexes2 = [(epsp_xindexes2[0] - 1)] + epsp_xindexes2
            t = [data['t'][x] - plot_start for x in xindexes]
            t2 = [data['t'][x] - data['t'][epsp_xindexes[0]] for x in epsp_xindexes]
            t3 = [data['t'][x] - data['t'][epsp_xindexes2[0]] for x in epsp_xindexes2]

            data_vmem_soma = [data['soma']['v'][x] for x in xindexes]
            data_vmem_syn = [data['synapse'][0]['v'][x] for x in xindexes]
            data_nr2a = [data['synapse'][0]['g_nr2a'][x]*1e3 for x in xindexes]
            data_nr2b = [data['synapse'][0]['g_nr2b'][x]*1e3 for x in xindexes]
            data_epsp1 = [data['soma']['v'][x] for x in epsp_xindexes]
            data_epsp2 = [data['soma']['v'][x] for x in epsp_xindexes2]
            
            epsp_threshold = 0.005
            for i in range(len(data_epsp1)-2):
                change1 = data_epsp1[i+1] - data_epsp1[i]
                change2 = data_epsp1[i+2] - data_epsp1[i+1]
                if (change2 - change1) > epsp_threshold:
                    epsp_start_time1 = t2[i+2]
                    break
            for i in range(len(data_epsp2)-2):
                change1 = data_epsp2[i+1] - data_epsp2[i]
                change2 = data_epsp2[i+2] - data_epsp2[i+1]
                if (change2 - change1) > epsp_threshold:
                    epsp_start_time2 = t3[i+2]
                    break
            t2 = [x - epsp_start_time1 + 5 for x in t2]
            t3 = [x - epsp_start_time2 + 5 for x in t3]

            hist_data = [x['weight'][-1] for x in data['synapse']]

            vmemid = 0
            nmdarid = 1
            weightid = 2
            epspid = 3


            N, bins, patches = axes[idx][weightid].hist(hist_data, color = 'gray', bins=np.arange(0.0, 2.6, 0.1), rwidth=0.9, label=r"${\omega}$", orientation='vertical')
            for i in range(len(N)):
                r = (2.6 - patches[i].get_x())/2.6
                g = patches[i].get_x()/2.6
                b = 0
                c = (r, g, b)
                patches[i].set_facecolor(c)

            axes[idx][vmemid].plot(t, data_vmem_syn, c = colors[0]/255, label = "$V_{d}$", linewidth=linewidth)
            axes[idx][vmemid].plot(t, data_vmem_soma, c = colors[1]/255, label = "$V_{s}$", linewidth=linewidth)
            axes[idx][nmdarid].plot(t, data_nr2a, c = colors[2]/255, label = "GluN2A", linewidth=linewidth)
            axes[idx][nmdarid].plot(t, data_nr2b, c = colors[3]/255, label = "GluN2B", linewidth=linewidth)
            axes[idx][epspid].plot(t2, data_epsp1, c = colors[0]/255, label = "$V_{before}$", linewidth=linewidth)
            axes[idx][epspid].plot(t3, data_epsp2, c = colors[1]/255, label = "$V_{after}$", linewidth=linewidth)

            axes[idx][weightid].set_xlim(0, 2.6)
            axes[idx][weightid].set_xticks([0, 1.0, 2.0, 2.6])
            axes[idx][vmemid].set_xlim(100, 1600)
            axes[idx][nmdarid].set_xlim(100, 1600)
            axes[idx][epspid].set_xlim(0, epsp_xlim)
            axes[idx][weightid].set_ylim(0, 50)
            axes[idx][vmemid].set_ylim(vmem_ylims[protocol_idx][0], vmem_ylims[protocol_idx][1])
            axes[idx][nmdarid].set_ylim(glu_ylims[protocol_idx][0], glu_ylims[protocol_idx][1])
            axes[idx][epspid].set_ylim(epsp_ylims[protocol_idx][0], epsp_ylims[protocol_idx][1])

            axes[idx][weightid].set_ylabel("Synapses #", fontsize=label_size)
            axes[idx][vmemid].set_ylabel("mV", fontsize=label_size)
            axes[idx][nmdarid].set_ylabel("nS", fontsize=label_size)
            axes[idx][epspid].set_ylabel("mV", fontsize=label_size)
            
            axes[idx][weightid].tick_params(labelsize=tick_size)
            axes[idx][vmemid].tick_params(labelsize=tick_size)
            axes[idx][nmdarid].tick_params(labelsize=tick_size)
            axes[idx][epspid].tick_params(labelsize=tick_size)

            axes[idx][vmemid].legend(fontsize=legend_size, loc=1)
            axes[idx][nmdarid].legend(fontsize=legend_size, loc=1)
            axes[idx][epspid].legend(fontsize=legend_size, loc=1)
            

            text_y = ((epsp_ylims[protocol_idx][1] - epsp_ylims[protocol_idx][0]) / 2) + epsp_ylims[protocol_idx][0]
            text_x = 30 + 3
            axes[idx][epspid].text(text_x, text_y, f"{labels[idx]}", fontsize=label_size2, ha='left', va='center').set_clip_on(False)

        axes[0][weightid].set_title("C) Weight\n", fontdict={'fontsize': title_size, 'fontweight': 'heavy'})
        axes[0][vmemid].set_title("A) Membrane potential\n", fontdict={'fontsize': title_size, 'fontweight': 'heavy'})
        axes[0][nmdarid].set_title("B) NMDAR conductance\n", fontdict={'fontsize': title_size, 'fontweight': 'heavy'})
        axes[0][epspid].set_title("D) EPSP\n", fontdict={'fontsize': title_size, 'fontweight': 'heavy'})
        axes[-1][weightid].set_xlabel("Weight", fontsize=label_size)
        axes[-1][vmemid].set_xlabel("ms", fontsize=label_size)
        axes[-1][nmdarid].set_xlabel("ms", fontsize=label_size)
        axes[-1][epspid].set_xlabel("ms", fontsize=label_size)

        filenamespl = filename.split('_')
        filenamespl.insert(-1, protocol_name)
        filenamespl = "_".join(filenamespl)

        savepath = f"{path}/{filenamespl}"
        print(f"Saving figure {savepath}")
        plt.savefig(savepath, dpi=dpi, bbox_inches='tight')
        if savemore:
            filenamespl = savemore.replace('vmem', f'vmem_{protocol_name}')
            filenamespl = filenamespl.split('.')
            filenamespl.insert(1, f"{protocol_idx+1}")
            filenamespl = ".".join(filenamespl)
            savepath = f"{path}/{filenamespl}"
            print(f"Saving figure {savepath}")
            plt.savefig(savepath, dpi=dpi, bbox_inches='tight')
        plt.close()


def plot_epsp_experiment(datas, filename, datanames, config, savemore=None, epspticks=None):
    fig, ax = plt.subplots(
        nrows=1, ncols=2,
        figsize=(24, 10),
        )

    labels = config['labels']['epsp']
    colors = config['colors']['epsp']
    colors = np.array(colors)

    title_size = 34
    tick_size = 22
    label_size = 32
    annotate_size = 20
    legend_size = 30
    
    titles = [
        "A) LTP",
        "B) LTD",
    ]
    for protocol in [0, 1]:
        weights = [datas[data][protocol]['epsp'] for data in datanames]
        pps = ax[protocol].bar(range(len(weights)), weights, tick_label=labels, color=colors/255)

        if epspticks:
            for p in pps:
                height = p.get_height()
                ax[protocol].annotate(f"{height:.2f}",
                    xy=(p.get_x() + p.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=annotate_size)

        ax[protocol].set_ylabel("EPSP (%)", fontsize=label_size)
        ax[protocol].set_title(titles[protocol], fontdict={'fontsize': title_size, 'fontweight': 'heavy'})

        ax[protocol].tick_params(labelsize=tick_size)
        ax[protocol].set_ylim(0, 220)
        ax[protocol].set_xticks(range(len(weights)), labels, rotation=30, ha='right')
    plt.tight_layout()
    savepath = f"{path}/{filename}"
    print(f"Saving figure {savepath}")
    plt.savefig(savepath, dpi=dpi)
    if savemore:
        savepath = f"{path}/{savemore}"
        print(f"Saving figure {savepath}")
        plt.savefig(savepath, dpi=dpi)
    plt.close()


def plot_epsp_nr2b(datas, filename, datanames, config, configname, savemore=None, epspticks=None, inset=False):
    fig, ax = plt.subplots(
        nrows=1, ncols=2,
        figsize=(34, 10),
        )

    titles = [
        "A) LTP",
        "B) LTD",
    ]
    
    labels = []
    nr2b_proportions = [config['runs'][x]['neuron_parameters']['nr2bmult'] for x in datanames]
    nr2b_aicd_multiplier = [config['runs'][x]['alzheimers'][1] for x in datanames]
    for lid in range(int(len(nr2b_proportions)/2)):
        nr2b = nr2b_proportions[int(len(nr2b_proportions)/2)+lid] * nr2b_aicd_multiplier[int(len(nr2b_proportions)/2)+lid]
        labels.append(f"{nr2b_proportions[lid]:.2f} | {nr2b:.2f}")
    title = config['titles'][configname]
    midpoint = int(len(datanames) / 2)

    colors = config['colors'][configname]
    colors = np.array(colors)


    title_size = 34
    tick_size = 22
    label_size = 32
    annotate_size = 20
    legend_size = 30

    
    for protocol in [0, 1]:
        control_weights = [datas[data][protocol]['epsp'] for data in datanames[:midpoint]]
        aicd_weights = [datas[data][protocol]['epsp'] for data in datanames[midpoint:]]

        width = 0.4
        pps = ax[protocol].bar((np.arange(len(labels)) + 1 - width).tolist(), control_weights, [width] * len(labels), color=colors[0]/255, label='Control')

        if epspticks:
            for p in pps:
                height = p.get_height()
                ax[protocol].annotate(f"{height:.2f}",
                    xy=(p.get_x() + p.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=annotate_size)
        pps = ax[protocol].bar((np.arange(len(labels)) + 1).tolist(), aicd_weights, [width] * len(labels), color=colors[1]/255, label=title)

        if epspticks:
            for p in pps:
                height = p.get_height()
                ax[protocol].annotate(f"{height:.2f}",
                    xy=(p.get_x() + p.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=annotate_size)

        ax[protocol].set_ylabel("EPSP (%)", fontsize=label_size)
        ax[protocol].set_xlabel("Proportion of active GluN2B-NMDAR", fontsize=label_size)

        ax[protocol].set_xticks((np.arange(len(labels)) + 1 - (width/2)).tolist(), labels)

        ax[protocol].set_title(titles[protocol], fontdict={'fontsize': title_size, 'fontweight': 'heavy'})
        ax[protocol].set_ylim(0, 220)

        ax[protocol].tick_params(labelsize=tick_size)
        ax[protocol].legend(fontsize=legend_size, loc=2, ncols=2)

    if inset:
        newax = ax[0].inset_axes([0.6, 0.6, 0.4, 0.4], anchor='NE')
        im = mpimg.imread('figures/elife.png')
        newax.imshow(im)
        newax.xaxis.set_tick_params(labelbottom=False)
        newax.yaxis.set_tick_params(labelleft=False)

        newax.set_xticks([])
        newax.set_yticks([])

    plt.tight_layout()
    savepath = f"{path}/{filename}"
    print(f"Saving figure {savepath}")
    plt.savefig(savepath, dpi=dpi)
    if savemore:
        savepath = f"{path}/{savemore}"
        print(f"Saving figure {savepath}")
        plt.savefig(savepath, dpi=dpi)
    plt.close()
