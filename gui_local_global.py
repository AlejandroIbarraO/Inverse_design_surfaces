import polyscope as ps
import polyscope.imgui as psim
import numpy as np
import tkinter as tk
import matplotlib.pyplot as plt
from tkinter import filedialog
from scipy import spatial
from igl import igl
from local_global import local_global
from phasor import phasor
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
from skimage import measure
root = tk.Tk()
root.withdraw()

lg = None
my_phasor = None
is_true1 = False
is_true2 = True
ui_int = 7
ui_float1 = -3.2
ui_float2 = 0.8
ui_color3 = (1., 0.5, 0.5)
ui_color4 = (0.3, 0.5, 0.5, 0.8)
ui_angle_rad = 0.2
ui_text = "some input text"
ui_options = ["Fabric", "Baromophs", "Direct Control"]
ui_options_selected = ui_options[0]
ui_lm1 = 1.25
ui_lm2 = 0.8944271909999159
ui_period = 2.0
def display_mesh():
    if lg is not None:
        ps_mesh = ps.register_surface_mesh("original", lg.v, lg.f, smooth_shade=True)
    
def lams_histogram(lams,ax= None,lams_target=None):
    x= np.linspace(0.5,2,200)
    if ax is None:
        _, ax = plt.subplots()
    ax.hist(lams[:,0], bins=x, density = True)
    ax.hist(lams[:,1], bins=x, density = True)
    if lams_target is not None:
        ax.vlines(lams_target,0,1)

    return ax
# Define our callback function, which Polyscope will repeatedly execute while running the UI.
# We can write any code we want here, but in particular it is an opportunity to create ImGui 
# interface elements and define a custom UI.
def callback():

    global is_true1, is_true2, ui_int, ui_float1, ui_float2, ui_color3, ui_color4, ui_text, ui_options_selected, ui_angle_rad

    global lg,ui_lm1,ui_lm2,ui_period,my_phasor

    psim.PushItemWidth(100)


    # == Show text in the UI

    psim.TextUnformatted("Simple Local-Global minimization parametrization")
    #psim.TextUnformatted("An important value: {}".format(42))
    #psim.Separator()


    # == Buttons

    if(psim.Button("Load OBJ")):        
        file_path = filedialog.askopenfilename(filetypes=[("3D mesh file", ".obj .stl")])
        lg = local_global(file_path)
        display_mesh()

    # By default, each element goes on a new line. Use this 
    # to put the next element on the _same_ line.
    psim.SameLine() 
    if(psim.Button("Precomputation")):
        if lg is not None:
            lg.precomputation()

    psim.PushItemWidth(200)
    changed = psim.BeginCombo("Pick one", ui_options_selected)
    if changed:
        for val in ui_options:
            _, selected = psim.Selectable(val, ui_options_selected==val)
            if selected:
                ui_options_selected = val
        psim.EndCombo()
    psim.PopItemWidth()

    changed, ui_lm1= psim.InputFloat("lambda 1", ui_lm1) 
    psim.SameLine() 
    changed, ui_lm2= psim.InputFloat("lambda 2", ui_lm2)

    changed, ui_int = psim.InputInt("N iterations", ui_int, step=1, step_fast=10) 
    if ui_int<0:
        ui_int = 0
    if(psim.Button("One Step L_G")):
        if lg is not None:

            if ui_options_selected == 'Direct control':
                lg.lam1 = ui_lm1
                lg.lam2 = ui_lm2
            else:
                lg.flavor = ui_options_selected
            lg.local_global_step()
    psim.SameLine() 
    if(psim.Button("N Step L_G")):
        _,lms,_ = lg.strech_and_angle()
        ax = lams_histogram(lms)
        if lg is not None:
            if ui_options_selected == 'Direct control':
                lg.lam1 = ui_lm1
                lg.lam2 = ui_lm2
            else:
                lg.flavor = ui_options_selected
            lg.local_global_n_step(ui_int)
        _,lms,_ = lg.strech_and_angle()
        ax = lams_histogram(lms,ax = ax,lams_target=[ui_lm1,ui_lm2])
        plt.show()
    if(psim.Button("Set fix point")):
        if lg is not None:
            mesh_name, id = ps.get_selection()
            lg.constrain(np.array([id]))
            point = lg.v[id].reshape(-1,3)
            fix_point = ps.register_point_cloud("fix point", point)

    if(psim.Button("Current uv map")):
        if lg is not None:
            uv_mesh = ps.register_surface_mesh("current uv", lg.u, lg.f, smooth_shade=True)
            centers,lms,angle_field = lg.strech_and_angle()
            index = is_true1*0+is_true2*1
            print(index)
            uv_mesh.add_scalar_quantity("angle", lms[:,index], defined_on='faces')
            uv_cloud = ps.register_point_cloud("centers", centers)
            vecs = np.array([np.cos(angle_field),np.sin(angle_field)]).T
            uv_cloud.add_vector_quantity("direction", vecs, enabled=True)
    
    psim.SameLine() 
    changed, is_true1 = psim.Checkbox("Lam1", is_true1) 
    if(changed): # optionally, use this conditional to take action on the new value
        if is_true1:    
            is_true2 = False
        else:
            is_true2 = True
    psim.SameLine() 
    changed, is_true2 = psim.Checkbox("Lam2", is_true2) 
    if(changed): # optionally, use this conditional to take action on the new value
        if is_true2:    
            is_true1 = False
        else:
            is_true1 = True
    changed, ui_period= psim.InputFloat("phasor period", ui_period) 
    if(psim.Button("Pattern")):
        if lg is not None:
            # create and array of phasors in a square latice inside the bounding box of the flatten configuration
            ulim = [lg.u[:,0].min(),lg.u[:,0].max()]
            vlim = [lg.u[:,1].min(),lg.u[:,1].max()]
            pu,pv = np.meshgrid(np.arange(ulim[0],ulim[1]+ui_period/2,ui_period/2),np.arange(vlim[0],vlim[1]+ui_period/2,ui_period/2))
            centers = np.array([pu.flatten(),pv.flatten(),np.zeros_like(pv.flatten())]).T
            l = igl.boundary_loop(lg.f)
            poly = Polygon(lg.u[l].tolist())
            index =[poly.contains(Point(cen[0],cen[1])) for cen in centers]
            centers = centers[index,:]
            uv_phasors = ps.register_point_cloud("centers-phasors", centers)
            tri_center,angle = lg.angles()
            tree = spatial.KDTree(tri_center)
            _,nieg = tree.query(centers[:,:2],1)
            

            vx = np.cos(angle)
            vy = np.sin(angle)
            vp_x = vx[nieg]
            vp_y = vy[nieg]
            vp = np.array([vp_x,vp_y]).T
            uv_phasors.add_vector_quantity("direction", vp, enabled=True)
            angle_phasor = np.arctan2(vp_y,vp_x)
            my_phasor = phasor()
            l = igl.avg_edge_length(lg.v, lg.f)*2.5
            #index = np.random.choice(np.arange(len(angle_face)),n_phasors)
            my_phasor.phase_sync = True
            if ui_options_selected == 'Fabric':
                my_phasor.phasor_compute(centers[:,:2],angle_phasor,1/ui_period,ui_period)
            else:
                my_phasor.phasor_compute(centers[:,:2],angle_phasor+np.pi/2,1/ui_period,ui_period)

            gabor =  np.sin(my_phasor.eval(lg.u,ui_period))
            uv_mesh = ps.register_surface_mesh("phasor", lg.u, lg.f, smooth_shade=True)
            #_,angle_field = lg.angles()
            uv_mesh.add_scalar_quantity("gabor_noise", gabor)
            subsampling = ui_period/4
            puq,pvq = np.meshgrid(np.arange(ulim[0]-subsampling,ulim[1]+2*subsampling,subsampling),np.arange(vlim[0]-subsampling,vlim[1]+2*subsampling,subsampling))
            centers = np.array([puq.flatten(),pvq.flatten()]).T
            image = np.zeros_like(puq).flatten()
            index =[poly.contains(Point(cen[0],cen[1])) for cen in centers]
            image[index] = np.sin(my_phasor.eval(centers[index],ui_period))*0.5+0.5
            image = (image-image.min())/(image.max()-image.min())
            image = image.reshape(puq.shape)
            fig, ax = plt.subplots()
            ax.imshow(image)
            contours = measure.find_contours(image, 0.5)
            np.save('data.npy', image)

            plt.show()

    
ps.init() 
ps.set_user_callback(callback)
ps.show()