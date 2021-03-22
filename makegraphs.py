from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt

last_index= 50;

live_count= np.zeros( last_index+1 );
apop_v_count= np.zeros( last_index+1 );
pyro_v_count= np.zeros( last_index+1 );
pyro_by_count= np.zeros( last_index+1 );

times = np.zeros( last_index+1 );

for n in range( 0,last_index+1 ):

	filename='output'+"%08i"%n+'.xml'

	#mcds=pyMCDS(filename,'output_a50_p50_Ril1b_300')
	mcds=pyMCDS(filename,'output/')
	times[n]= mcds.get_time()

	# Flags

	cycle=mcds.data['discrete_cells']['total_volume']
	pyro_flag=mcds.data['discrete_cells']['cell_pyroptosis_flag']
	pyro_by_flag=mcds.data['discrete_cells']['cell_bystander_pyroptosis_flag']
	apo_virus_flag=mcds.data['discrete_cells']['cell_virus_induced_apoptosis_flag']
	pyro_time=mcds.data['discrete_cells']['cell_pyroptosis_time']
	apo_time=mcds.data['discrete_cells']['cell_apo_time']
	
	# Cell is viable:
	live = np.argwhere((cycle==2494) & (pyro_flag==0) & (apo_virus_flag==0) & (pyro_by_flag==0)).flatten()

	# Cell is apoptosing due to viral load:
	apop_v = np.argwhere((apo_virus_flag==1)  & (apo_time==0)).flatten()


	# Cell is pyroptosing due to viral load
	pyro_v = np.argwhere( (pyro_flag==1) & (pyro_time==1) ).flatten()
	# Cell is bystander pyroptosing
	pyro_by = np.argwhere( (pyro_by_flag==1)& (pyro_time==1) ).flatten()

	live_count[n] = len(live)
	apop_v_count[n] = apop_v_count[n-1]+len(apop_v)
	pyro_v_count[n] = pyro_v_count[n-1]+len(pyro_v)
	pyro_by_count[n] = pyro_by_count[n-1]+len(pyro_by)

plt.plot(live_count, color='blue', label='viable')
plt.plot(apop_v_count, color='black', label='apoptosing (viral load)')
plt.plot(pyro_v_count, color='orange', label='pyroptosing (viral load)')
plt.plot(pyro_by_count, color='red', label='pyroptosing (bystander)')


plt.ylabel('Number of cells', size=20 )

plt.xlabel('Time (hrs)', size=20 )

plt.legend(loc= 'upper right')

plt.show()