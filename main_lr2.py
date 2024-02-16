import re

import numpy as np


def parse_liberty(filename, gate_name):
  # Local variables
  fl_str_in_trans = 0
  fl_col_out_cap = 0
  fl_gate_start = 0
  fl_gate_cap = 0
  fl_gate_fall = 0
  fl_gate_rise = 0
  fl_gate_fall_trans = 0
  fl_gate_rise_trans = 0

  # Table headers
  in_trans_list = []
  out_cap_list = []

  # Gate capacitance
  gate_cap = 0

  # Gate analysis
  fall_tab = []
  rise_tab = []
  fall_trans_tab = []
  rise_trans_tab = []

  # Parse file
  with open(filename, 'r', encoding='utf-8') as infile:
    for line in infile:
      # Fetch transition time
      if (fl_str_in_trans):
        in_trans_list = np.array(line.replace("\n", "").split(","),
                                 dtype='float')
        fl_str_in_trans = 0

      # Fetch output capacitance
      if (fl_col_out_cap):
        out_cap_list = np.array(line.replace("\n", "").split(","),
                                dtype='float')
        fl_col_out_cap = 0

      # Fetch gate capacitance
      if (fl_gate_cap):
        gate_cap = float(line)
        fl_gate_cap = 0

      # Fetch gate fall
      if (fl_gate_fall == 1):
        if (line != '\n'):
          fall_tab.append(
              np.array(line.replace("\n", "").split(","), dtype='float'))
        else:
          fl_gate_fall += 1

      # Fetch gate rise
      if (fl_gate_rise == 1):
        if (line != '\n'):
          rise_tab.append(
              np.array(line.replace("\n", "").split(","), dtype='float'))
        else:
          fl_gate_rise += 1

      # Fetch gate fall transition
      if (fl_gate_fall_trans == 1):
        if (line != '\n'):
          fall_trans_tab.append(
              np.array(line.replace("\n", "").split(","), dtype='float'))
        else:
          fl_gate_fall_trans += 1

      # Fetch gate rise transition
      if (fl_gate_rise_trans == 1):
        if (line != '\n'):
          rise_trans_tab.append(
              np.array(line.replace("\n", "").split(","), dtype='float'))
        else:
          fl_gate_rise_trans += 1

      # Check gate EOL
      if (fl_gate_start and fl_gate_fall == 2 and fl_gate_rise == 2
          and fl_gate_fall_trans == 2 and fl_gate_rise_trans == 2):
        fl_gate_start = 0
        fl_gate_cap = 0
        fl_gate_fall = 0
        fl_gate_rise = 0
        fl_gate_fall_trans = 0
        fl_gate_rise_trans = 0

      # Check input transition
      if (re.match('string_input_transtion', line)):
        fl_str_in_trans = 1
      # Check output capacitance
      if (re.match('column_output_capacitance', line)):
        fl_col_out_cap = 1
      # Check gate tables start
      if (re.match(gate_name, line)):
        fl_gate_start = 1
      # Check gate capacitance
      if (fl_gate_start and re.match('capacitance', line)):
        fl_gate_cap = 1
      # Check gate fall
      if (fl_gate_start and re.match('cell_fall', line)):
        fl_gate_fall = 1
      # Check gate rise
      if (fl_gate_start and re.match('cell_rise', line)):
        fl_gate_rise = 1
      # Check gate fall transition
      if (fl_gate_start and re.match('fall_transition', line)):
        fl_gate_fall_trans = 1
      # Check gate rise transition
      if (fl_gate_start and re.match('rise_transition', line)):
        fl_gate_rise_trans = 1
  return in_trans_list, out_cap_list, gate_cap, fall_tab, rise_tab, fall_trans_tab, rise_trans_tab


def find_borders(in_data, in_data_list, in_closest_index):
  # Get left and right index
  if (in_data >= in_data_list[in_closest_index]):
    ind_l = in_closest_index
    ind_r = in_closest_index + 1
  else:
    ind_l = in_closest_index - 1
    ind_r = in_closest_index
  return ind_l, ind_r


def find_borders_value(i_l, i_r, list):
  return list[i_l], list[i_r]


def bilinear_interpol(table, str_l, str_r, col_l, col_r, hdr_str_l, hdr_str_r,
                      hdr_col_l, hdr_col_r, in_str, in_col):
  aprox_str_col_l = ((hdr_col_r - in_col) /
                     (hdr_col_r - hdr_col_l) * table[str_l][col_l] +
                     (in_col - hdr_col_l) /
                     (hdr_col_r - hdr_col_l) * table[str_l][col_r])
  aprox_str_col_r = ((hdr_col_r - in_col) /
                     (hdr_col_r - hdr_col_l) * table[str_r][col_l] +
                     (in_col - hdr_col_l) /
                     (hdr_col_r - hdr_col_l) * table[str_r][col_r])
  aprox_result = ((hdr_str_r - in_str) /
                  (hdr_str_r - hdr_str_l) * aprox_str_col_l +
                  (in_str - hdr_str_l) /
                  (hdr_str_r - hdr_str_l) * aprox_str_col_r)
  return aprox_result


class Gate:
  """Cell transition and delay class"""
  cells_count = 0

  def __init__(self, lib_file, gate_name):
    self.name = gate_name
    self.trans_header = []
    self.cap_header = []
    self.capacitance = 0
    self.fall_tab = []
    self.rise_tab = []
    self.fall_trans_tab = []
    self.rise_trans_tab = []
    (self.trans_header, self.cap_header, self.capacitance, self.fall_tab,
     self.rise_tab, self.fall_trans_tab,
     self.rise_trans_tab) = parse_liberty(lib_file, gate_name)
    Gate.cells_count += 1

  def get_delay(self, edge, in_tran, out_cap):
    # Check out of range input
    if ((in_tran < self.trans_header[0] or in_tran > self.trans_header[-1])
        or (out_cap < self.cap_header[0] or out_cap > self.cap_header[-1])
        or not ((edge == "p") or (edge == "n"))):
      print("Out of range")
    else:
      # Get closest table index
      in_trans_ind = np.abs(self.trans_header - in_tran).argmin()
      out_cap_ind = np.abs(self.cap_header - out_cap).argmin()
      # Find table border
      tr_i_l, tr_i_r = find_borders(in_tran, self.trans_header, in_trans_ind)
      cap_i_l, cap_i_r = find_borders(out_cap, self.cap_header, out_cap_ind)
      tr_v_l, tr_v_r = find_borders_value(tr_i_l, tr_i_r, self.trans_header)
      cap_v_l, cap_v_r = find_borders_value(cap_i_l, cap_i_r, self.cap_header)
      # Interpolate value
      if (edge == "p"):
        linear = self.rise_tab[in_trans_ind][out_cap_ind]
        bilinear = bilinear_interpol(self.rise_tab, tr_i_l, tr_i_r, cap_i_l,
                                     cap_i_r, tr_v_l, tr_v_r, cap_v_l, cap_v_r,
                                     in_tran, out_cap)
      else:
        linear = self.fall_tab[in_trans_ind][out_cap_ind]
        bilinear = bilinear_interpol(self.fall_tab, tr_i_l, tr_i_r, cap_i_l,
                                     cap_i_r, tr_v_l, tr_v_r, cap_v_l, cap_v_r,
                                     in_tran, out_cap)
      return linear, bilinear

  def get_transition(self, edge, in_tran, out_cap):
    # Check out of range input
    if ((in_tran < self.trans_header[0] or in_tran > self.trans_header[-1])
        or (out_cap < self.cap_header[0] or out_cap > self.cap_header[-1])
        or not ((edge == "p") or (edge == "n"))):
      print("Out of range")
    else:
      # Get closest table index
      in_trans_ind = np.abs(self.trans_header - in_tran).argmin()
      out_cap_ind = np.abs(self.cap_header - out_cap).argmin()
      # Find table border
      tr_i_l, tr_i_r = find_borders(in_tran, self.trans_header, in_trans_ind)
      cap_i_l, cap_i_r = find_borders(out_cap, self.cap_header, out_cap_ind)
      tr_v_l, tr_v_r = find_borders_value(tr_i_l, tr_i_r, self.trans_header)
      cap_v_l, cap_v_r = find_borders_value(cap_i_l, cap_i_r, self.cap_header)
      # Interpolate value
      if (edge == "p"):
        linear = self.rise_trans_tab[in_trans_ind][out_cap_ind]
        bilinear = bilinear_interpol(self.rise_trans_tab, tr_i_l, tr_i_r,
                                     cap_i_l, cap_i_r, tr_v_l, tr_v_r, cap_v_l,
                                     cap_v_r, in_tran, out_cap)
      else:
        linear = self.fall_trans_tab[in_trans_ind][out_cap_ind]
        bilinear = bilinear_interpol(self.fall_trans_tab, tr_i_l, tr_i_r,
                                     cap_i_l, cap_i_r, tr_v_l, tr_v_r, cap_v_l,
                                     cap_v_r, in_tran, out_cap)
      return linear, bilinear

  def get_propagate(self, edge, in_tran, out_cap):
    del_lin, del_bilin = self.get_delay(edge, in_tran, out_cap)
    tran_lin, tran_bilin = self.get_transition(edge, in_tran, out_cap)
    return del_lin, del_bilin, tran_lin, tran_bilin


# Input data
gate_trans_type = "p"
gate_in_trans = "16"
gate_out_cap = "43"

# Convert input
gate_in_trans = float(gate_in_trans)
gate_out_cap = float(gate_out_cap)

# Add objects
g_and = Gate("lib1.lib", "AND")
g_xnor = Gate("lib1.lib", "XNOR")
g_nor = Gate("lib1.lib", "NOR")
print(Gate.cells_count)

#Solve task
task_delay = 0
task_tran = 0
_, tmp_del, _, tmp_tran = g_and.get_propagate("p", gate_in_trans,
                                              g_xnor.capacitance)
task_delay += tmp_del
_, tmp_del, _, tmp_tran = g_xnor.get_propagate("n", tmp_tran,
                                               g_nor.capacitance)
task_delay += tmp_del
_, tmp_del, _, tmp_tran = g_nor.get_propagate("p", tmp_tran, gate_out_cap)
task_delay += tmp_del
print("Total delay:", task_delay)
print("Transition:", tmp_tran)