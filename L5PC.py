import os
import sys
import neuron
from neuron import gui
from neuron import h
import cell_parent as cp
import models
import utils
from matplotlib import pyplot as plt
import pickle as pkl
import copy

class L5PC_model(cp.cell_parent,):
    def __init__(self, branch_num = 78, nseg = 50, dend_L = 200, Hay = False):
        super().__init__()
        self.dendrites = []
        h.load_file("nrngui.hoc")
        
        if not Hay:
            h.load_file("import3d.hoc")
            h("objref cell, tobj") #neuron object    
            morphology_file = "morphologies/cell1.asc"
            h.load_file('allen_model.hoc')
            h.execute("cell = new allen_model()")  # replace?
            nl = h.Import3d_Neurolucida3()
            nl.quiet = 1
            nl.input(morphology_file)
            i3d = h.Import3d_GUI(nl, 0)
            i3d.instantiate(h.cell)
            self.hoc_cell=h.cell
            self.hoc_cell.geom_nseg()
            self.hoc_cell.biophys()
            dend = self.hoc_cell.dend[branch_num]
    #         self.hoc_cell.delete_axon()
        if Hay:
            morphology_file = "morphologies/cell1.asc"
            h.load_file("nrngui.hoc")
            h("objref cvode")
            h("cvode = new CVode()")
            h("cvode.active(0)")
            h.load_file("import3d.hoc") 
            h.load_file("models/L5PCbiophys3.hoc")
            h.load_file("models/L5PCtemplate.hoc")
            self.hoc_cell = h.L5PCtemplate(morphology_file)
            dend = self.hoc_cell.dend[branch_num]
            dend.L = dend_L
            dend.nseg = nseg
            dend.diam = 0.75
            dend.Ra = 150
            self.dend = dend 
            
        self.soma = self.hoc_cell.soma[0]
        
#         dend.L = dend_L
#         dend.nseg = nseg
#         dend.diam = 1
#         self.dend = dend 
        self.dendrites.append(dend) #comment back in !!!
        self.setValidLocations()
        self.save_file_name = 'mech inits/mechinit branch num_'+str(branch_num)+'_ nseg_'+str(nseg)+'_.pkl'
        self.model_name = 'L5PC'
        if Hay:
            if os.path.exists(self.save_file_name):
                self.fih = neuron.h.FInitializeHandler(1, self.load_saved_state)
            else:
                self.run_sim_with_save()
        
        self.heads = []
        self.necks = []
        
        'to add branches to Hay model, comment in following lines'
#         for sec_num in range(1,15):
#             dend = h.Section(name="dend_sec_"+str(sec_num))
#             dend.L = dend_L#should maybe be random?
# #             dend.diam = dend_diam
#             dend.Ra = 100
#             dend.cm = 2 #1?
#             dend.nseg = dend_L//5
#             dend.insert('pas')
# #             dend.insert('Ih')
# #             dend.gIhbar_Ih = 0.0002
#             dend.g_pas =  0.0000467
# #             dend.e_pas = -90
#             self.dendrites.append(dend)
#              
#         for sec_num in range(1,len(self.dendrites)):
#             if sec_num == 0:
#                 self.dendrites[0].connect(self.soma, 1, 0)
#             else:
#                 connect_to = utils.find_dend_sec_path(sec_num)[-2]
#                 self.dendrites[sec_num].connect(self.dendrites[connect_to], 1, 0)
#                 
#         print(self.dendrites)
        self.get_dend_TRs(), self.get_soma_TRs()
        
    def run_sim(self, rec_locations, sim_time, dt):
        super().run_sim(rec_locations, sim_time, dt,v_init=-70)

    def run_sim_with_save(self, rec_locations = [0.1], sim_time = 400, dt = 0.1):
        self.run_sim(rec_locations, sim_time, dt)
        self.save_state()
        plt.plot(self.t, self.soma_v)
        plt.show()
        
#     def save_state(self, file_name):
#         f = h.File(file_name)
#         svstate = h.SaveState()
#         svstate.save()
#         svstate.fwrite(f)
#         f.close()

    def setValidLocations(self):
        self.trees = []
        self.trees.append(self.hoc_cell.apic)
        self.trees.append(self.hoc_cell.dend)
        self.trees.append(self.hoc_cell.soma)
        self.trees.append(self.hoc_cell.axon)
        self.allLocations = []
        for i in range(4):
            treeNum = i
            tree = self.trees[treeNum]
            for j in range(len(tree)):
                sec = tree[j]
                for k in range(sec.nseg):
                    self.allLocations.append((treeNum,j,k))

    def loc2Seg(self, loc):
        tree = self.trees[loc[0]]
        section = tree[loc[1]]
#         seg = section(self.getSegLoc(loc))
        seg = section(loc[2])
        return seg
    
    def get_seg_from_thruple(self, key):
        tree = self.trees[key[0]]        
        section = tree[key[1]]
        nsegs = float(section.nseg)
        segNum = key[2]
        segLoc = (segNum / nsegs + 0.5 / nsegs)
        return section(segLoc)
    
    def save_state(self):
        mechLocMap = {}
        for loc in self.allLocations:
            seg_num = self.get_seg_from_thruple(loc)
            seg = self.get_seg_from_thruple(loc)
            mechMap = {}
            mechMap['v'] = seg.v
            for mech in seg:
                if '__' not in mech.name():
                    for var in dir(mech):
                        if '__' not in var and 'name' not in var and 'next' not in var and 'is_ion' not in var and 'segment' not in var:
                            try:
                                tm = eval('seg.' + var)
                                mechMap[var] = (None, tm)
                                #print('try', var, tm, type(tm))
                            except:
                                tm = eval('seg.' + str(mech) + '.'+ var)
                                mechMap[var] = (str(mech), tm)
                                #print('except', var, tm, type(tm))
            mechLocMap[loc] = mechMap.copy()
            #print(list(mechLocMap))
        fn = self.save_file_name
        with open(fn, 'wb') as f:
            pkl.dump(mechLocMap, f)
        
    def load_saved_state(self):
#         print('a')
#         f = h.File()
#         f.close()
#         print('b')
#         f.ropen(file_name)
#         print('b1')
#         svstate = h.SaveState()
#         print('c')
#         f.seek(0)
#         print('c1')
#         svstate.fread(f)
#         print('d')
#         h.finitialize()
#         print('e')
#         svstate.restore()
#         t = 0 # t is one of the "states"
# #         if (cvode.active()) {
# #             cvode.re_init()
# #           } else {
# #             fcurrent()
# #           }
# #         frecord_init()
#         print('f')
#         svstate.restore()
#         print('g')

#         self.data = pkl.load(open(file_name, "rb" ))
#         for key in self.data.keys():
#             seg = self.get_seg_from_thruple(key)
#             #seg.v = data[key]['v']
#             for key_1 in self.data[key].keys():
#                 print(str(key_1))
#                 print('seg.'+str(key_1)+' = '+str(self.data[key][key_1]))
#                 exec('seg.'+str(key_1)+' = '+str(self.data[key][key_1]))
        mechLocMap = pkl.load(open(self.save_file_name, "rb" ))
        arr = [mechLocMap[loc]['v'] for loc in self.allLocations if loc[0] == 2]
        for loc in mechLocMap.keys():
            seg = self.get_seg_from_thruple(loc)
            mechMap = mechLocMap[loc]
            seg.sec.v = mechMap['v']
            seg.v = mechMap['v']
            for mech in seg:
                if '__' not in mech.name():
                    for var in dir(mech):
                        if '__' not in var and 'name' not in var and 'next' not in var and 'is_ion' not in var and 'segment' not in var:
                            tup = mechMap[var]
                            if tup[0] == None:
                                exec('seg.' + str(var) + ' = ' + str(tup[1]))
                            else:
                                exec('seg.' + str(tup[0]) + '.' + str(var) + ' = ' + str(tup[1]))
            
    def loc2Seg(self, loc):
        treeNum = loc[0]
        tree = self.trees[treeNum]
        section = tree[loc[1]]
        seg = section(self.getSegLoc(loc))
        return seg
     
if __name__ == "__main__": 
# cell = L5PC_model(file_name = 'savedmechs.pkl')
    cell = L5PC_model(Hay = True)
    cell.plot_TRs()
# #Assertion failed: i < prl_->count(), file ../../../nrn/src/nrniv/../nrncvode/netcvode.cpp, line 5619
#     cell.plot_TRs()
# #plt.show()
# #cell.save_state('test')
# #cell.run_sim_with_save([0.1], 5000, 0.025, 'steady state')
# cell.run_sim([0.1], 50, 0.025,)
#cell.load_saved_state()
