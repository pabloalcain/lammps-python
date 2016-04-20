import matplotlib
matplotlib.use('Agg')
params = {'axes.labelsize': 12,
          'font.size': 15,
          'legend.fontsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15,
          'text.usetex': False,
          'lines.linewidth': 2.0,
          'lines.markersize': 10.0,
          'legend.numpoints': 1,
          'axes.color_cycle': ('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3',
                               '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3')}

#COLOR = ('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e',
#         '#e6ab02', '#a6761d', '#666666')


def _format(axis):
  """
  Format axis to make them pretty
  """
  axis.spines['top'].set_visible(False)
  axis.spines['right'].set_visible(False)
  axis.spines['left'].set_visible(False)
  axis.get_xaxis().tick_bottom()
  axis.get_yaxis().tick_left()
  axis.tick_params(axis='x', direction='out')
  axis.tick_params(axis='y', length=0)
  for spine in axis.spines.values():
    spine.set_position(('outward', 5))
  axis.grid(axis='y', color='0.9', linestyle='-', linewidth=1)
  axis.set_axisbelow(True)

matplotlib.rcParams.update(params)
