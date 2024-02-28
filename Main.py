import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
import random
from tqdm import tqdm

from Constants import *
from Extra import *
from LineManager import *
from Physics_Tranfers import *
from Satellite import *

from panda3d.core import AmbientLight , DirectionalLight , Point3 , MouseButton
from panda3d.core import Vec3 , KeyboardButton , TextureStage , TransparencyAttrib
from panda3d.core import LightAttrib , NodePath , CardMaker , NodePath , TextNode
from panda3d.core import AntialiasAttrib, loadPrcFileData

from direct.showbase.ShowBase import ShowBase
from direct.task import Task
from direct.gui.OnscreenText import OnscreenText
from direct.gui.OnscreenImage import OnscreenImage




class MyApp(ShowBase):

    def __init__(self):
        ShowBase.__init__(self)
        self.anti_antialiasing(is_on=True)

        # OTV
        otv_init_elements = orbital_elements(inclination = 20,
                                             raan        = 0,
                                             eccentricity= 0.01,
                                             arg_perigee = 0,
                                             mean_anomaly= 60,
                                             a           = 2*R_earth)
        self.otv , self.otv_node = self.make_object(elements=otv_init_elements)
        
        # Target
        target_init_elements = orbital_elements(inclination = 10,
                                                raan        = 200,
                                                eccentricity= 0.01,
                                                arg_perigee = 0,
                                                mean_anomaly= 50,
                                                a           = 3*R_earth)
        self.target , self.target_node = self.make_object(elements=target_init_elements)


        # Transfers:
        # Inclination
        self.do_transfer_inc = True
        self.transfer_d_inc = self.otv.elements.inclination - self.target.elements.inclination
        self.show_transfer_inc_line = False
        self.inc_boost_to_log = 0

        # RAAN
        self.do_transfer_raan = True
        self.raan_boost_to_log = 0

        # Hohmann
        self.do_hohmann = True
        self.show_hohmann_lines = False
        self.boost_to_log = 0


        self.setup_inc_transfer_params()
        self.setup_raan_transfer_params()
        self.setup_hohmann_transfer_params()


        self.setup_scene()
        self.taskMgr.add(self.check_keys, "check_keys_task")
        self.accept("space" , self.on_space_pressed)
        self.game_is_paused = False
        self.accept("a" , self.on_a_pressed) # Show planes
        # To check the function: elements -> pos,vel -> elements
        self.otv_true , self.otv_true_node = self.make_object(elements=otv_init_elements)


        self.visualise = False

        if self.visualise:
            self.taskMgr.doMethodLater(1/DT, self.renderer, 'renderer')
        else:
            for _ in tqdm(range(N_log)):
                self.renderer(None)


        

        


    def setup_scene(self):
        self.setup_log_arrays()
        self.setup_skybox()
        self.setup_camera()
        self.setup_nodes()
        self.setup_lights()
        self.setup_hud()


    def renderer(self, task):

        if not self.game_is_paused:
            self.otv.update(dt=DT)
            self.target.update(dt=DT)

            
            
            rotate_object(self.earth , [0.05 , 0 , 0])
            rotate_object(self.cloud , [0.05 , 0 , 0])
            self.env_visual_update()
            
            
            if self.do_transfer_inc:
                self.make_inc_transfer()
            if self.do_transfer_raan and self.inc_boost_done:
                self.make_raan_transfer()
            if self.do_hohmann and self.raan_boost_done:
                self.make_hohmann_transfer()

            self.log_values()

            if abs(self.otv.elements.mean_anomaly - 360) < 2:
                self.on_space_pressed()
                

        return Task.cont
        

    

    def log_values(self):
        if self.can_log == False:
            return
        
        self.log_timer -= 1
        if self.log_timer <= 0:
            self.log_timer = self.log_interval
        else:
            return
        

        # otv:
        self.log_otv_rad.append(np.linalg.norm(self.otv.position / R_earth))
        self.log_otv_inc.append(self.otv.elements.inclination)
        self.log_otv_raan.append(self.otv.elements.raan)
        self.log_otv_nu.append(self.otv.elements.mean_anomaly_360)

        if self.boost_to_log > 0:
            self.boost_to_log = 0
            self.log_hoh_boost.append(np.linalg.norm(self.otv.position / R_earth))
        else:
            self.log_hoh_boost.append(None)

        if self.inc_boost_to_log > 0:
            self.inc_boost_to_log = 0
            self.log_inc_boost.append(self.otv.elements.inclination)
        else:
            self.log_inc_boost.append(None)

        if self.raan_boost_to_log > 0:
            self.raan_boost_to_log = 0
            self.log_raan_boost.append(self.otv.elements.raan)
        else:
            self.log_raan_boost.append(None)

        self.log_otv_ang_vel.append(np.degrees(self.otv.find_angular_velocity()))

        _ , u_final , _ = algorithm_40(otv=self.otv , target=self.target)
        nu_boost = np.degrees(u_final) - self.otv.elements.arg_perigee
        nu_otv = self.otv.elements.mean_anomaly_360
        diff = abs(nu_boost - nu_otv)
        
        self.log_raan_u.append(min(diff, 360 - diff))

        # target:
        self.log_target_rad.append(np.linalg.norm(self.target.position / R_earth))
        self.log_target_inc.append(self.target.elements.inclination)
        self.log_target_raan.append(self.target.elements.raan)
        self.log_target_nu.append(self.target.elements.mean_anomaly_360)
        self.log_target_ang_vel.append(np.degrees(self.target.find_angular_velocity()))

        self.log_otv_target_dist.append(np.linalg.norm(self.otv.position - self.target.position) / R_earth)



        if len(self.log_otv_rad) >= N_log:
            self.can_log = False
            self.plot_values()


    def setup_log_arrays(self):

        # otv:
        self.log_otv_rad = []
        self.log_otv_inc = []
        self.log_otv_raan = []
        self.log_otv_nu = []
        self.log_hoh_boost = []
        self.log_inc_boost = []
        self.log_raan_boost = []
        self.log_otv_ang_vel = []

        self.log_raan_u = []

        # target:
        self.log_target_rad = []
        self.log_target_inc = []
        self.log_target_raan = []
        self.log_target_nu = []
        self.log_target_ang_vel = []

        self.log_otv_target_dist = []



        self.log_interval = 1 # frames between each log
        self.log_timer = self.log_interval
        self.can_log = True
        

    def plot_values(self):
        x = range(len(self.log_otv_rad))
        fig, axs = plt.subplots(2, 3, figsize=(15, 5))

        # Radius of orbit plot
        axs[0,0].plot(x, self.log_otv_rad, 'r' , label='otv')
        axs[0,0].plot(x, self.log_target_rad, 'b--' , label='target')
        axs[0,0].scatter(x, self.log_hoh_boost, color='purple' , label='boost')
        axs[0,0].set_title('Radius of orbit')
        axs[0,0].set_ylabel('R_earth')
        axs[0,0].legend()
        axs[0,0].set_ylim(0, 4)

        # Inclination plot
        axs[0,1].plot(x, self.log_otv_inc, 'r' , label='otv')
        axs[0,1].plot(x, self.log_target_inc, 'b--' , label='target')
        axs[0,1].scatter(x, self.log_inc_boost, color='purple' , label='boost')
        axs[0,1].set_title('Inclination')
        axs[0,1].set_ylabel('Degrees')
        axs[0,1].legend()

        # RAAN plot
        axs[0,2].plot(x, self.log_otv_raan, 'r' , label='otv')
        axs[0,2].plot(x, self.log_target_raan, 'b--' , label='target')
        axs[0,2].scatter(x, self.log_raan_boost, color='purple' , label='boost')
        axs[0,2].set_title('Raan')
        axs[0,2].set_ylabel('Degrees')
        axs[0,2].legend()

        # OTV -> Target distance plot
        axs[1,0].plot(x , self.log_otv_target_dist , 'r')
        axs[1,0].set_title('OTV Target distance')
        axs[1,0].set_ylabel('R_earth')
        axs[1,0].set_ylim(0 , 5)

        # True anomaly plot
        axs[1,1].plot(x , self.log_otv_nu, 'r' , label='otv')
        axs[1,1].plot(x , self.log_target_nu, 'b--' , label='target')
        axs[1,1].set_title("True anomaly")
        axs[1,1].set_ylabel('Degrees')
        axs[1,1].legend()

        # Angular velocity plot
        # axs[1,2].plot(x , self.log_otv_ang_vel , 'r' , label='otv')
        # axs[1,2].plot(x , self.log_target_ang_vel , 'b--' , label='target')
        # axs[1,2].set_title("Angular velocity")
        # axs[1,2].set_ylabel('Degrees / ?')
        axs[1,2].plot(x , self.log_raan_u)
        axs[1,2].legend()

        plt.tight_layout()
        plt.show()

        

    def setup_raan_transfer_params(self):
        self.initial_delay_raan = 0 # To replace with the phase time to the correct anomaly
        self.raan_boost_done = False
        

    def make_raan_transfer(self):
        # Renderer function
        if self.raan_boost_done:
            return
        

        dv_raan_only , u_final , theta = algorithm_40(otv=self.otv , target=self.target)
        nu_boost = np.degrees(u_final) - self.otv.elements.arg_perigee # should check the range of arg_perigee
        thr = 1.0
        nu_otv = self.otv.elements.mean_anomaly

        diff = abs(nu_boost - nu_otv)
        
        if min(diff, 360 - diff) < thr:
            pass
        else:
            return
        
        if self.raan_boost_done == False:
            self.setup_hohmann_transfer_params() # inc -> raan -> hohmann

            self.raan_boost_done = True
            self.raan_boost_to_log = 1

            print(f"Do raan boost, otv true anomaly={self.otv.elements.mean_anomaly_360}")
            
            rad_dir = normalize_vector(self.otv.position)
            vel_dir = normalize_vector(self.otv.velocity)
            top_dir = np.cross(vel_dir , rad_dir)

            rotation = R.from_rotvec(rad_dir * (-theta/2))
            boost_dir = rotation.apply(top_dir)

            self.otv.velocity += boost_dir * dv_raan_only


    

    def setup_inc_transfer_params(self):
        phase_to_0 = simple_phase(object=self.otv , target_anomaly=0)
        #phase_to_180 = simple_phase(object=self.otv , target_anomaly=180)
        
        self.initial_delay_inc = phase_to_0 # min(phase_to_0 , phase_to_180)
        print(f"Inc: t_0={phase_to_0}")# , t_180={phase_to_180}")

        self.inc_boost_done = False


    def make_inc_transfer(self):
        # Renderer function
        if self.inc_boost_done:
            return

        # if self.initial_delay_inc > 0:
        #     self.initial_delay_inc -= DT
        #     if self.initial_delay_inc < 0:
        #         self.initial_delay_inc = 0
        #     return

        nu = self.otv.elements.mean_anomaly_360
        thr = 1
        # if nu < thr or (360-nu)< thr:
        #     print(nu)
        # else:
        #     return
        
        if abs(180-nu) < thr:
            pass#print(nu)
        else:
            return
        
        if self.inc_boost_done == False:
            self.setup_raan_transfer_params() # inc -> raan -> hohmann

            self.inc_boost_done = True
            self.inc_boost_to_log = 1
            print(f"Do inc boost, otv true anomaly={self.otv.elements.mean_anomaly_360} , {self.otv.elements.arg_perigee}")
            

            ang_ = -np.deg2rad(self.transfer_d_inc/2)

            rad_dir = normalize_vector(self.otv.position)

            # Create a rotation object from the axis and angle
            rotation = R.from_rotvec(rad_dir * ang_)

            vel_dir = normalize_vector(self.otv.velocity)
            top_dir = np.cross(vel_dir , rad_dir)

            # Apply rotation to the vector
            boost_dir = rotation.apply(top_dir)


            inc_1 = self.otv.elements.inclination
            inc_2 = inc_1 + self.transfer_d_inc

            self.otv.velocity += inc_dv(v1=np.linalg.norm(self.otv.velocity) , inc1=inc_1 , inc2=inc_2) * boost_dir




    def setup_hohmann_transfer_params(self):
        self.hohmann_a1 = self.otv.elements.a
        self.hohmann_a2 = self.target.elements.a
        #print(self.hohmann_a1/R_earth , self.hohmann_a2/R_earth)
        self.hoh_dv1 , self.hoh_dv2 = hohmann_dv(self.hohmann_a1 , self.hohmann_a2) # a in [m]
        self.hoh_dt = hohmann_time(self.hohmann_a1 , self.hohmann_a2)

        #print(f"Hohmann: dv1={self.hoh_dv1} , dv2={self.hoh_dv2} , t={self.hoh_dt}")

        self.initial_delay = algorithm_45(otv=self.otv , target=self.target)
        self.boost_1_done = False
        self.boost_2_done = False

        

        
    def make_hohmann_transfer(self):
        # Renderer function
        if self.boost_1_done and self.boost_2_done:
            return

        if self.initial_delay > 0:
            self.initial_delay -= DT
            if self.initial_delay < 0:
                self.initial_delay = 0
            return

        if self.boost_1_done == False:
            print(f"Hohmann dv1, at r={np.linalg.norm(self.otv.position)/R_earth}")
            self.boost_1_done = True
            self.boost_to_log = 1

            v_dir = normalize_vector(self.otv.velocity)
            dv = self.hoh_dv1 * v_dir
            self.otv.velocity += dv
            return
        
        # Boost 1 done
        
        if self.hoh_dt > 0:
            self.hoh_dt -= DT
            return
        
        # Hohmann time done
        
        if self.boost_2_done == False and self.boost_1_done == True:
            print(f"Hohmann dv2, at r={np.linalg.norm(self.otv.position)/R_earth}")
            self.boost_2_done = True
            self.boost_to_log = 1

            v_dir = normalize_vector(self.otv.velocity)
            dv = self.hoh_dv2 * v_dir
            self.otv.velocity += dv
            return
        
        # Boost 2 done

        
        

    def env_visual_update(self):
        otv_pos = self.otv.position / R_earth # Normalised position
        self.otv_node.setPos(otv_pos[0] , otv_pos[1] , otv_pos[2])

        target_pos = self.target.position / R_earth
        self.target_node.setPos(target_pos[0] , target_pos[1] , target_pos[2])


        self.otv_true.elements = self.otv.elements
        true_pos = np.array(self.otv_true.find_pos_vel()[0]) / R_earth
        self.otv_true_node.setPos(true_pos[0] , true_pos[1] , true_pos[2])
        #print(f"pos diff: {np.linalg.norm(otv_pos - true_pos)}")

        self.update_otv_trail()
        self.update_hud()


        line_pos = []
        for pos in self.otv_trail_nodes:
            if pos != (0,0,0):
                line_pos.append(pos)

        if self.show_hohmann_lines == False:
            self.line_manager.update_line('trail_line', line_pos, color=(1, 0, 0, 1))

        vel_value = np.linalg.norm(self.otv.velocity)

        vel_point_1 = tuple(otv_pos)
        vel_point_2 = tuple(vel_point_1 + 0.5*self.otv.velocity/vel_value)
        self.line_manager.update_line('vel_line', [vel_point_1, vel_point_2], color=(0, 1, 0, 1))

        if self.show_hohmann_lines == False:
            rad_point_1 = tuple(otv_pos + normalize_vector(-otv_pos)*0.5)
            self.line_manager.update_line("radial_line" , [vel_point_1 , rad_point_1] , color=(0,0,1,1))

        if self.show_hohmann_lines:
            centre_point = tuple([0,0,0])

            # t0
            otv_point = tuple(normalize_vector(otv_pos) * 3.5)
            target_point = tuple(normalize_vector(self.target.position / R_earth) * 3.5)

            self.line_manager.update_line('hohmann_otv' , [centre_point , otv_point] , color=(0, 1, 0, 1))
            self.line_manager.update_line('hohmann_target' , [centre_point , target_point] , color=(0, 1, 0, 1))

            
        
        if self.show_transfer_inc_line:
            # d_inc = self.transfer_inc2 - self.transfer_inc1
            ang_ = -np.deg2rad(self.transfer_d_inc/2)

            rad_dir = normalize_vector(self.otv.position)

            # Create a rotation object from the axis and angle
            rotation = R.from_rotvec(rad_dir * ang_)

            vel_dir = normalize_vector(self.otv.velocity)
            top_dir = np.cross(vel_dir , rad_dir)

            # Apply rotation to the vector
            boost_dir = rotation.apply(top_dir)

            inc_v_point_2 = tuple(self.otv.position + 10*boost_dir/vel_value)
            self.line_manager.update_line("inc_vel_line" , [vel_point_1 , inc_v_point_2] , color=(1,0,1,1))

            new_v_inc_point_2 = tuple(vel_point_2+2*boost_dir/vel_value)
            self.line_manager.update_line("new_vel_inc" , [vel_point_1 , new_v_inc_point_2] , color=(1,1,1,1))


    def setup_nodes(self):
        self.make_earth()

        self.all_debris_nodes = []
        self.all_debris_planets = []

        self.otv_trail_nodes = []
        self.otv_trail_counter = 10
        self.make_otv_trail()

        self.line_manager = LineManager(self.render)
        self.line_manager.make_line('trail_line', [(0, 0, 0), (0, 0, 0)], color=(1, 0, 0, 1))
        self.line_manager.make_line('vel_line', [(0, 0, 0), (0, 0, 0)], color=(1, 0, 0, 1))
        self.line_manager.make_line('radial_line', [(0, 0, 0), (0, 0, 0)], color=(0, 0, 1, 1))
        #if self.show_transfer_inc_line:
        self.line_manager.make_line('inc_vel_line' , [(0, 0, 0), (0, 0, 0)], color=(0, 0, 1, 1))
        self.line_manager.make_line('new_vel_inc' , [(0, 0, 0), (0, 0, 0)], color=(0, 0, 1, 1))

        self.line_manager.make_line('hohmann_otv' , [(0, 0, 0), (0, 0, 0)], color=(0, 0, 1, 1))
        self.line_manager.make_line('hohmann_target' , [(0, 0, 0), (0, 0, 0)], color=(0, 0, 1, 1))
        self.line_manager.make_line('hohmann_phase_ang' , [(0, 0, 0), (0, 0, 0)], color=(0, 0, 1, 1))

        self.setup_planes_visualisation()
        




    def make_object(self , elements):
        init_pos = [1,1,1]
        node = self.make_sphere(size=0.03 , low_poly=True)
        node.setPos(init_pos[0] , init_pos[1] , init_pos[2])
        node.reparentTo(self.render)

        sat = Satellite(elements)
        pos = sat.position
        node.setPos(pos[0] , pos[1] , pos[2])

        return sat , node
    


    def make_otv_trail(self):
        init_pos = (0,0,0) # At first, all the nodes are hidden inside the earth
        n = 100

        for _ in range(n):
            node = init_pos
            self.otv_trail_nodes.append(node)


    def update_otv_trail(self):
        # update every frame
        self.otv_trail_counter -= 1

        if self.otv_trail_counter <= 0:
            self.otv_trail_counter = 4
        else:
            return
        
        self.otv_trail_nodes[0] = self.otv_node.getPos()
        self.otv_trail_nodes = self.otv_trail_nodes[1:] + [self.otv_trail_nodes[0]]


        
    def make_earth(self):
        self.earth = self.make_sphere(size=1)
        self.earth.setPos(0, 0, 0)
        self.earth.setHpr(0, 90, 0)

        # Load textures
        base_texture = self.loader.loadTexture("textures/earth2.jpg")
        gloss_texture = self.loader.loadTexture("textures/earth_spe.png")
        glow_texture = self.loader.loadTexture("textures/earth_glow.png")

        # Base Texture Stage
        base_ts = TextureStage('base')
        base_ts.setMode(TextureStage.MModulate)
        self.earth.setTexture(base_ts, base_texture)

        # Gloss Texture Stage
        gloss_ts = TextureStage('gloss')
        gloss_ts.setMode(TextureStage.MGloss)
        self.earth.setTexture(gloss_ts, gloss_texture)

        # Glow Texture Stage
        glow_ts = TextureStage('glow')
        glow_ts.setMode(TextureStage.MAdd)
        self.earth.setTexture(glow_ts, glow_texture)

        self.earth.reparentTo(self.render)


        self.cloud = self.make_sphere(1.01)
        self.cloud.setPos(0, 0, 0)
        self.cloud.setHpr(0 , 90 , 0)
        cloud_texture = self.loader.loadTexture("textures/earth_cloud.png")

        #cloud_texture.setFormat(Texture.F_srgb_alpha)
        self.cloud.setTexture(cloud_texture)

        # Enable transparency for the object
        self.cloud.setTransparency(TransparencyAttrib.M_alpha)

        # You may need to adjust the transparency value based on your texture
        self.cloud.setAlphaScale(0)
        self.cloud.reparentTo(self.render)
        


        

    def make_a_debris(self):
        # Planet object
        elements = orbital_elements(inclination=random.uniform(0 , 90),
                                        raan=30,
                                        eccentricity=random.uniform(0.0 , 0.3),
                                        arg_perigee=0,
                                        mean_anomaly=0,
                                        a=random.uniform(1.8 , 2.2))
        
        sat_object = Satellite(elements)
        self.all_debris_planets.append(sat_object)

        # Node
        position = sat_object.position
        node = self.make_sphere(0.01 , low_poly=True)
        node.setPos(position[0] , position[1] , position[2])
        node.reparentTo(self.render)
        self.all_debris_nodes.append(node)

        


    def setup_lights(self):
        # Ambient light
        ambient_light = AmbientLight("ambient_light")
        ambient_light.setColor((0.1, 0.1, 0.1, 1))
        ambient_light_node = self.render.attachNewNode(ambient_light)
        self.render.setLight(ambient_light_node)

        # Directional light
        directional_light = self.render.attachNewNode("directional_light")
        d_light = DirectionalLight("d_light")
        intensity = 4
        d_light.setColor((intensity, intensity, intensity, 1))
        d_light_np = directional_light.attachNewNode(d_light)
        d_light_np.setHpr(45, -15, 0)
        self.render.setLight(d_light_np)


    def setup_camera(self):
        self.rotation_speed = 50.0
        self.elevation_speed = -50.0

        self.distance_to_origin = 10.0
        self.distance_speed = 0.1
        self.min_dist = 4
        self.max_dist = 16

        self.angle_around_origin = 0.0
        self.elevation_angle = 0.0
        self.last_mouse_x = 0
        self.last_mouse_y = 0
        
        self.camera.setPos(0, -10, 0)
        self.camera.lookAt(0, 0, 0)

        self.camLens.setNear(0.1)
        self.camLens.setFar(100.0)
        
        self.disableMouse() # Enable mouse control for the camera
        self.accept("mouse1", self.mouse_click)

        self.taskMgr.add(self.update_camera_task, "update_camera_task")

    def setup_planes_visualisation(self):

        # Axis planes
        s = (4,4)
        c = (0 , 0.51 , 0.71 , 0.3)
        self.visualisation_plane_1 = self.create_plane(size=s , position=(0,0,0) , rotation=(90,0,0) , color=c)
        self.visualisation_plane_2 = self.create_plane(size=s , position=(0,0,0) , rotation=(0,90,0) , color=c)
        self.visualisation_plane_3 = self.create_plane(size=s , position=(0,0,0) , rotation=(0,0,90) , color=c)

        self.visualisation_plane_1.hide()
        self.visualisation_plane_2.hide()
        self.visualisation_plane_3.hide()

        self.visualisation_plane_1.setAttrib(LightAttrib.makeAllOff())
        self.visualisation_plane_2.setAttrib(LightAttrib.makeAllOff())
        self.visualisation_plane_3.setAttrib(LightAttrib.makeAllOff())

        self.visualisation_plane_is_on = False

        

    def toggle_plane_visualisation(self):
        self.visualisation_plane_is_on = not self.visualisation_plane_is_on

        if self.visualisation_plane_is_on:
            self.visualisation_plane_1.show()
            self.visualisation_plane_2.show()
            self.visualisation_plane_3.show()
        else:
            self.visualisation_plane_1.hide()
            self.visualisation_plane_2.hide()
            self.visualisation_plane_3.hide()
    

    def create_plane(self, size, position, rotation, color=(1, 1, 1, 1)):
        card_maker = CardMaker('plane')
        w, h = size
        card_maker.setFrame(-w / 2, w / 2, -h / 2, h / 2)  # Set the size of the plane
        plane_np = NodePath(card_maker.generate())
        plane_np.reparentTo(self.render)
        plane_np.setPos(position)  # Set position
        plane_np.setHpr(rotation)  # Set rotation

        # Set color and alpha
        r, g, b, a = color  # Unpack the color tuple
        plane_np.setColor(r, g, b, a)  # Apply color and alpha
        plane_np.setTwoSided(True)
        plane_np.setTransparency(TransparencyAttrib.MAlpha)

        return plane_np
        

    def setup_skybox(self):
        size = 20
        distance = 20

        texture_list = []
        for i in range(1,7):
            texture_list.append(self.loader.loadTexture(f"skybox/Skybox_{i}.jpg"))
        
        plane_1 = self.loader.loadModel("models/plane.obj")
        plane_2 = self.loader.loadModel("models/plane.obj")
        plane_3 = self.loader.loadModel("models/plane.obj")
        plane_4 = self.loader.loadModel("models/plane.obj")
        plane_5 = self.loader.loadModel("models/plane.obj")
        plane_6 = self.loader.loadModel("models/plane.obj")

        base_ts = TextureStage('base_ts')
        base_ts.setMode(TextureStage.MReplace)

        plane_1.setTexture(base_ts , texture_list[0])
        plane_1.setScale(size)
        plane_2.setTexture(base_ts , texture_list[0])
        plane_2.setScale(size)
        plane_3.setTexture(base_ts , texture_list[0])
        plane_3.setScale(size)
        plane_4.setTexture(base_ts , texture_list[0])
        plane_4.setScale(size)
        plane_5.setTexture(base_ts , texture_list[0])
        plane_5.setScale(size)
        plane_6.setTexture(base_ts , texture_list[0])
        plane_6.setScale(size)

        plane_1.setPos(0, 0, distance)
        plane_1.setHpr(0 , -90 , 0)

        plane_2.setPos(0, 0, -distance)
        plane_2.setHpr(0 , 90 , 0)

        plane_3.setPos(distance, 0, 0)
        plane_3.setHpr(90 , 0 , 0)

        plane_4.setPos(-distance, 0, 0)
        plane_4.setHpr(-90 , 0 , 0)

        plane_5.setPos(0, distance, 0)
        plane_5.setHpr(-180 , 0 , 0)

        plane_6.setPos(0, -distance, 0)
        plane_6.setHpr(0 , 0 , 90)

        plane_1.reparentTo(self.render)
        plane_2.reparentTo(self.render)
        plane_3.reparentTo(self.render)
        plane_4.reparentTo(self.render)
        plane_5.reparentTo(self.render)
        plane_6.reparentTo(self.render)


    



    def setup_hud(self):
        y_st = 0.9
        y_sp = 0.1
        x_po = -1.3

        # Top left orbital params
        self.title_1 = self.add_text_label(text="Orbital Elements:" , pos=(x_po , y_st) , scale=0.08)
        self.label_1 = self.add_text_label(text="label 1" , pos=(x_po , y_st-y_sp))
        self.label_2 = self.add_text_label(text="label 2" , pos=(x_po , y_st-y_sp*2))
        self.label_3 = self.add_text_label(text="label 3" , pos=(x_po , y_st-y_sp*3))
        self.label_4 = self.add_text_label(text="label 4" , pos=(x_po , y_st-y_sp*4))
        self.label_5 = self.add_text_label(text="label 5" , pos=(x_po , y_st-y_sp*5))
        self.label_6 = self.add_text_label(text="label 6" , pos=(x_po , y_st-y_sp*6))

        # Top right pause
        self.pause_label = self.add_text_label(text="II" , pos=(0 , y_st))
        self.pause_label.hide()
        
        # Middle Right
        self.title_3 = self.add_text_label(text="Inclination change:" , pos=(-x_po , y_st) , scale=0.08 , alignment_mode=TextNode.ARight)
        self.inc_label_1 = self.add_text_label(text="" , pos=(-x_po , y_st-y_sp) , alignment_mode=TextNode.ARight)
        self.inc_label_2 = self.add_text_label(text="" , pos=(-x_po , y_st-y_sp*2) , alignment_mode=TextNode.ARight)
        self.inc_label_3 = self.add_text_label(text="" , pos=(-x_po , y_st-y_sp*3) , alignment_mode=TextNode.ARight)
        
        self.title_4 = self.add_text_label(text="Raan change:" , pos=(-x_po , y_st-y_sp*6) , scale=0.08 , alignment_mode=TextNode.ARight)
        self.raan_label_1 = self.add_text_label(text="" , pos=(-x_po , y_st-y_sp*7) , alignment_mode=TextNode.ARight)
        self.raan_label_2 = self.add_text_label(text="" , pos=(-x_po , y_st-y_sp*8) , alignment_mode=TextNode.ARight)
        self.raan_label_3 = self.add_text_label(text="" , pos=(-x_po , y_st-y_sp*9) , alignment_mode=TextNode.ARight)

        # Middle Left
        self.title_2 = self.add_text_label(text="Hohmann phasing:" , pos=(x_po , y_st-y_sp*9) , scale=0.08)
        self.label_8 = self.add_text_label(text="label 8" , pos=(x_po , y_st-y_sp*10))
        self.label_9 = self.add_text_label(text="label 9" , pos=(x_po , y_st-y_sp*11))
        self.label_10 = self.add_text_label(text="label 10" , pos=(x_po , y_st-y_sp*12))
        self.label_11 = self.add_text_label(text="label 11" , pos=(x_po , y_st-y_sp*13))


        y_st = -0.9
        # Bottom left constants
        self.label_G = self.add_text_label(text=f"G = {G} m³ km⁻¹ s⁻²" , pos=(x_po , y_st))
        self.label_M = self.add_text_label(text=f"M = {M} kg" , pos=(x_po , y_st+y_sp))
        self.label_J2 = self.add_text_label(text=f"J2 = {J2}" , pos=(x_po , y_st+y_sp*2))
        
        # Bottom right
        self.frames_left_label = self.add_text_label(text="_" , pos=(-x_po , y_st) , alignment_mode=TextNode.ARight)


    def update_hud(self):
        i = self.otv.elements.inclination
        raan = self.otv.elements.raan
        e = self.otv.elements.eccentricity
        arg_perigee = self.otv.elements.arg_perigee
        mean_anomaly = self.otv.elements.mean_anomaly_360
        a = self.otv.elements.a


        i = round(float(i) , 3)
        raan = round(float(raan))
        e = round(float(e) , 3)
        arg_perigee = int(arg_perigee)
        mean_anomaly = int(mean_anomaly)
        a = round(float(a)/R_earth , 3)
        

        self.label_1.setText(f"i = {i}°")
        self.label_2.setText(f"raan = {raan}°")
        self.label_3.setText(f"e = {e}")
        self.label_4.setText(f"arg_perigee = {arg_perigee}°")
        self.label_5.setText(f"true anomaly = {mean_anomaly}°")
        self.label_6.setText(f"a = {a} Rₑ")

        
        # Inc and Raan labels
        phase_to_0 = simple_phase(object=self.otv , target_anomaly=360)
        phase_to_180 = simple_phase(object=self.otv , target_anomaly=180)

        self.inc_label_1.setText(f"T min = {round(float(self.initial_delay_inc),4)}")
        self.inc_label_2.setText(f"T to 0 = {round(float(phase_to_0),2)}")
        self.inc_label_3.setText(f"T to 180 = {round(float(phase_to_180),2)}")

        phase_to_90 = simple_phase(object=self.otv , target_anomaly=90)
        phase_to_270 = simple_phase(object=self.otv , target_anomaly=270)

        self.raan_label_1.setText(f"T min = {round(self.initial_delay_raan , 4)}")
        self.raan_label_2.setText(f"T to 90 = {round(float(phase_to_90),2)}")
        self.raan_label_3.setText(f"T to 270 = {round(float(phase_to_270),2)}")

        # Algorithm 45 parts
        lead_angle , phase_angle , initial_phase_angle , T_wait = algorithm_45(otv=self.otv , target=self.target , prints=True , debug_msg=False)
        self.label_8.setText(f"Lead angle = {round(np.degrees(lead_angle))}°")
        self.label_9.setText(f"Phase angle = {round(np.degrees(phase_angle))}°")
        self.label_10.setText(f"Angle diff = {round(np.degrees(initial_phase_angle))}°")

        if self.inc_boost_done == False and self.do_transfer_inc:
            self.label_11.setText(f"T wait = not set")
        elif self.do_transfer_inc or not self.do_transfer_inc:
            self.label_11.setText(f"T wait = {round(self.initial_delay,2)}")
        
        
        pc_done = int(100*len(self.log_otv_rad)/N_log)
        self.frames_left_label.setText(f"{pc_done}%")
    


    
    def mouse_click(self):
        # Check if the left mouse button is down
        if self.mouseWatcherNode.isButtonDown(MouseButton.one()):
            # Enable mouse motion task
            self.last_mouse_x = self.mouseWatcherNode.getMouseX()
            self.last_mouse_y = self.mouseWatcherNode.getMouseY()
            self.taskMgr.add(self.update_camera_task, "update_camera_task")

    def update_camera_task(self, task):
        # Check if the left mouse button is still down
        if self.mouseWatcherNode.isButtonDown(MouseButton.one()):
            # Get the mouse position
            current_mouse_x = self.mouseWatcherNode.getMouseX()
            current_mouse_y = self.mouseWatcherNode.getMouseY()

            # Check if the mouse has moved horizontally
            if current_mouse_x != self.last_mouse_x:
                # Adjust the camera rotation based on the mouse horizontal movement
                self.angle_around_origin -= (current_mouse_x - self.last_mouse_x) * self.rotation_speed

            # Check if the mouse has moved vertically
            if current_mouse_y != self.last_mouse_y:
                # Adjust the camera elevation based on the mouse vertical movement
                self.elevation_angle += (current_mouse_y - self.last_mouse_y) * self.elevation_speed
                self.elevation_angle = max(-90, min(90, self.elevation_angle))  # Clamp the elevation angle

            self.update_camera_position()

            self.last_mouse_x = current_mouse_x
            self.last_mouse_y = current_mouse_y

            return task.cont
        else:
            # Disable the mouse motion task when the left button is released
            self.taskMgr.remove("update_camera_task")
            return task.done
        

    def update_camera_position(self):
        radian_angle = np.radians(self.angle_around_origin)
        radian_elevation = np.radians(self.elevation_angle)
        x_pos = self.distance_to_origin * np.sin(radian_angle) * np.cos(radian_elevation)
        y_pos = -self.distance_to_origin * np.cos(radian_angle) * np.cos(radian_elevation)
        z_pos = self.distance_to_origin * np.sin(radian_elevation)

        self.camera.setPos(Vec3(x_pos, y_pos, z_pos))
        self.camera.lookAt(Point3(0, 0, 0))
        

    def check_keys(self, task):
        # Check if the key is down
        if self.mouseWatcherNode.is_button_down(KeyboardButton.up()):
            self.move_forward()
        if self.mouseWatcherNode.is_button_down(KeyboardButton.down()):
            self.move_backward()
        
        return task.cont

    def move_forward(self):
        if self.distance_to_origin > self.min_dist:
            self.distance_to_origin -= self.distance_speed
            self.update_camera_position()


    def move_backward(self):
        if self.distance_to_origin < self.max_dist:
            self.distance_to_origin += self.distance_speed
            self.update_camera_position()
        

    def make_sphere(self , size=1 , low_poly=False):
        path = "models/sphere5.obj"
        if low_poly:
            path = "models/low_poly_sphere.obj"
        sphere = self.loader.loadModel(path)
        sphere.setScale(size)
        return sphere
    

    def add_text_label(self , text="PlaceHolder" , pos=(-1 , 1) , scale=0.06 , alignment_mode=TextNode.ALeft):
        custom_font = self.loader.loadFont('textures/SF-Pro.ttf')
        text_label = OnscreenText(text=text,
                                    pos=pos, # Position on the screen
                                    scale=scale, # Text scale
                                    fg=(1, 1, 1, 1), # Text color (R, G, B, A)
                                    bg=(0, 0, 0, 0.5), # Background color (R, G, B, A)
                                    align=alignment_mode, # Text alignment
                                    font=custom_font,
                                    mayChange=True) # Allow text to change dynamically
        return text_label
        
    def add_image(self, image_path, pos=(0, 0), scale=1, parent=None):
        if parent is None:
            parent = self.render2d  # Use self.render2d if no parent is specified.
        img = OnscreenImage(image=image_path, pos=pos, scale=scale, parent=parent)
        img.setTransparency(TransparencyAttrib.MAlpha)


    def on_a_pressed(self):
        self.toggle_plane_visualisation()

    def on_z_pressed(self):
        if self.visualisation_plane_4.isHidden():
            self.visualisation_plane_4.show()
        else:
            self.visualisation_plane_4.hide()

    def on_space_pressed(self):
        self.game_is_paused = not self.game_is_paused

        if self.game_is_paused:
            self.pause_label.show()
        else:
            self.pause_label.hide()

    def anti_antialiasing(self , is_on):
        if is_on:
            loadPrcFileData('', 'multisamples 4')  # Enable MSAA
            self.render.setAntialias(AntialiasAttrib.MAuto)
    


if __name__ == "__main__":
    app = MyApp()
    app.run()