import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
import random

from Constants import *
from Extra import *
from LineManager import *
from Physics_Tranfers import *
from Satellite import *

from panda3d.core import AmbientLight , DirectionalLight , Point3 , MouseButton
from panda3d.core import Vec3 , KeyboardButton , TextureStage , TransparencyAttrib
from panda3d.core import LightAttrib , NodePath , CardMaker , NodePath , TextNode

from direct.showbase.ShowBase import ShowBase
from direct.task import Task
from direct.gui.OnscreenText import OnscreenText
from direct.gui.OnscreenImage import OnscreenImage





class MyApp(ShowBase):

    def __init__(self):
        ShowBase.__init__(self)
        self.setup_scene()


        otv_init_elements = orbital_elements(inclination=0,
                                             raan=0,
                                             eccentricity=0,
                                             arg_perigee=0,
                                             mean_anomaly=0,
                                             a=2)
        self.otv , self.otv_node = self.make_object(elements=otv_init_elements)
        

        target_init_elements = orbital_elements(inclination=0,
                                                raan=0,
                                                eccentricity=0,
                                                arg_perigee=0,
                                                mean_anomaly=0,
                                                a=3)
        self.target , self.target_node = self.make_object(elements=target_init_elements)


        self.do_hohmann = True
        self.hohmann_a1 = 2
        self.hohmann_a2 = 3

        self.do_transfer_inc = False
        self.transfer_inc1 = 0
        self.transfer_inc2 = 10

        self.do_transfer_raan = False
        self.transfer_raan1 = 0
        self.transfer_raan2 = 10


        self.setup_hohmann_transfer_params()
        self.setup_inc_transfer_params()
        self.setup_raan_transfer_params()


        self.taskMgr.add(self.check_keys, "check_keys_task")
        self.taskMgr.doMethodLater(DT, self.renderer, 'renderer')


        self.accept("space" , self.on_space_pressed)
        self.game_is_paused = False
        self.accept("a" , self.on_a_pressed) # Show planes
        

        


    def setup_scene(self):
        self.setup_log_arrays()
        self.setup_skybox()
        self.setup_camera()
        self.setup_nodes()
        self.setup_lights()
        self.setup_hud()


    def renderer(self, task):

        if not self.game_is_paused:
            rotate_object(self.earth , [0.05 , 0 , 0])
            rotate_object(self.cloud , [0.05 , 0 , 0])

            self.otv.update(dt=DT)
            self.target.update(dt=DT)
            self.env_visual_update()
            
            if self.do_hohmann:
                self.make_hohmann_transfer()
            if self.do_transfer_inc:
                self.make_inc_transfer()
            if self.do_transfer_raan:
                self.make_raan_transfer()

            self.log_values()

        return Task.cont
        

    

    def log_values(self):
        if self.can_log == False:
            return
        
        self.log_timer -= 1
        if self.log_timer <= 0:
            self.log_timer = self.log_interval
        else:
            return
        
        
        i = self.otv.elements.inclination
        raan = self.otv.elements.raan


        e = self.otv.elements.eccentricity
        a = self.otv.elements.a
        inst_r = np.linalg.norm(self.otv.position)


        self.log_a.append(inst_r)
        self.log_a2.append(np.linalg.norm(self.target.position))
        self.log_e.append(e)
        self.log_i.append(i)
        self.log_raan.append(raan)


        if self.boost_to_log > 0:
            self.boost_to_log = 0
            self.log_hoh_boost.append(inst_r)
        else:
            self.log_hoh_boost.append(None)


        if len(self.log_a) >= N_log:
            self.can_log = False
            self.plot_values()


    def setup_log_arrays(self):
        self.log_a = []
        self.log_a2 = []
        self.log_e = []
        self.log_i = []
        self.log_raan = []

        self.log_hoh_boost = []

        self.log_interval = 1 # frames between each log
        self.log_timer = self.log_interval
        self.can_log = True
        

    def plot_values(self):

        x = range(len(self.log_a))
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))

        axs[0].plot(x, self.log_a, 'r' , label='a')
        #axs[0].plot(x , self.log_a2 , 'g' , label='a2')
        axs[0].scatter(x, self.log_hoh_boost, color='purple' , label='boost')
        axs[0].set_title('Semi-major axis')
        axs[0].legend()
        axs[0].set_ylim(0, 5)

        axs[1].plot(x, self.log_i, 'g' , label='inc')
        axs[1].set_title('Inclination')
        axs[1].legend()

        lg_raan = self.log_raan
        axs[2].plot(x, lg_raan, 'b' , label='raan')
        axs[2].set_title('Raan')
        axs[2].legend()

        plt.tight_layout()
        plt.show()

        

    
        

    
    

    def setup_raan_transfer_params(self):
        self.initial_delay_raan = 3
        self.raan_boost_done = False
        

    def make_raan_transfer(self):
        # Renderer function
        if self.raan_boost_done:
            return

        if self.initial_delay_raan > 0:
            self.initial_delay_raan -= DT
            return
        
        ang_thr = 5
        if abs(self.otv.elements.mean_anomaly-90) < ang_thr or abs(self.otv.elements.mean_anomaly-270) < ang_thr:
            self.raan_boost_done = True

            vel_dir = normalize_vector(self.otv.velocity)
            rad_dir = normalize_vector(-self.otv.position)
            boost_dir = np.cross(vel_dir , rad_dir)

            self.otv.velocity += raan_dv(v1=np.linalg.norm(self.otv.velocity) , raan1=self.transfer_raan1 , raan2=self.transfer_raan2) * boost_dir



    

    def setup_inc_transfer_params(self):
        self.initial_delay_inc = 4
        self.inc_boost_done = False

    def make_inc_transfer(self):
        # Renderer function
        if self.inc_boost_done:
            return

        if self.initial_delay_inc > 0:
            self.initial_delay_inc -= DT
            return
        
        ang_thr = 1
        if self.otv.elements.mean_anomaly < ang_thr or abs(self.otv.elements.mean_anomaly-180) < ang_thr:
            self.inc_boost_done = True

            inc_1 = self.transfer_inc1
            inc_2 = self.transfer_inc2
            d_inc = inc_2 - inc_1
            ang_ = -np.deg2rad(d_inc/2)

            rad_dir = normalize_vector(self.otv.position)

            # Create a rotation object from the axis and angle
            rotation = R.from_rotvec(rad_dir * ang_)

            vel_dir = normalize_vector(self.otv.velocity)
            top_dir = np.cross(vel_dir , rad_dir)

            # Apply rotation to the vector
            boost_dir = rotation.apply(top_dir)



            
            
            self.otv.velocity += inc_dv(v1=np.linalg.norm(self.otv.velocity) , inc1=inc_1 , inc2=inc_2) * boost_dir




    def setup_hohmann_transfer_params(self):
        self.hoh_dv1 , self.hoh_dv2 = hohmann_dv(self.hohmann_a1 , self.hohmann_a2)
        self.hoh_dt = hohmann_time(self.hohmann_a1 , self.hohmann_a2)

        self.initial_delay = phase_time(otv=self.otv , target=self.target) - np.pi/4
        self.boost_1_done = False
        self.boost_2_done = False

        self.boost_to_log = 0

        
    def make_hohmann_transfer(self):
        # Renderer function
        if self.boost_1_done and self.boost_2_done:
            return

        if self.initial_delay > 0:
            self.initial_delay -= DT
            return

        if self.boost_1_done == False:
            self.boost_1_done = True
            self.boost_to_log = 1
            print(self.otv.elements.a)

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
            self.boost_2_done = True
            self.boost_to_log = 1
            print(self.otv.elements.a)

            v_dir = normalize_vector(self.otv.velocity)
            dv = self.hoh_dv2 * v_dir
            self.otv.velocity += dv
            return
        
        # Boost 2 done

        
        

    def env_visual_update(self):

        otv_pos = self.otv.position
        self.otv_node.setPos(otv_pos[0] , otv_pos[1] , otv_pos[2])

        target_pos = self.target.position
        self.target_node.setPos(target_pos[0] , target_pos[1] , target_pos[2])

        self.update_otv_trail()
        self.update_hud()


        line_pos = []
        for pos in self.otv_trail_nodes:
            if pos != (0,0,0):
                line_pos.append(pos)


        self.line_manager.update_line('trail_line', line_pos, color=(1, 0, 0, 1))

        vel_value = np.linalg.norm(self.otv.velocity)

        vel_point_1 = tuple(self.otv.position)
        vel_point_2 = tuple(vel_point_1 + 0.5*self.otv.velocity/vel_value)
        self.line_manager.update_line('vel_line', [vel_point_1, vel_point_2], color=(0, 1, 0, 1))

        rad_point_1 = tuple(self.otv.position + normalize_vector(-self.otv.position)*0.5)
        self.line_manager.update_line("radial_line" , [vel_point_1 , rad_point_1] , color=(0,0,1,1))


        

        d_inc = self.transfer_inc2 - self.transfer_inc1
        ang_ = -np.deg2rad(d_inc/2)

        rad_dir = normalize_vector(self.otv.position)

        # Create a rotation object from the axis and angle
        rotation = R.from_rotvec(rad_dir * ang_)

        vel_dir = normalize_vector(self.otv.velocity)
        top_dir = np.cross(vel_dir , rad_dir)

        # Apply rotation to the vector
        boost_dir = rotation.apply(top_dir)

        inc_v_point_2 = tuple(self.otv.position + 0.5*boost_dir/vel_value)
        self.line_manager.update_line("inc_vel_line" , [vel_point_1 , inc_v_point_2] , color=(1,0,1,1))

        new_v_inc_point_2 = tuple(vel_point_2+0.5*boost_dir/vel_value)
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
        self.line_manager.make_line('inc_vel_line' , [(0, 0, 0), (0, 0, 0)], color=(0, 0, 1, 1))
        self.line_manager.make_line('new_vel_inc' , [(0, 0, 0), (0, 0, 0)], color=(0, 0, 1, 1))

        self.setup_planes_visualisation()
        




    def make_object(self , elements):
        init_pos = [1,1,1]
        node = self.make_sphere(size=0.03 , low_poly=True)
        node.setPos(init_pos[0] , init_pos[1] , init_pos[2])
        node.reparentTo(self.render)

        return Satellite(elements) , node
    


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


        # Orbital plane
        # self.visualisation_plane_4 = self.create_plane(size=s , position=(0,0,0) , rotation=(0,0,0) , color=(1 , 1 , 1 , 0.7))
        # self.update_orbital_plane()
        # self.visualisation_plane_4.setAttrib(LightAttrib.makeAllOff())
        # self.visualisation_plane_4.hide()

        
    # def update_orbital_plane(self):
    #     data = self.otv_planet.get_orbital_elements()

    #     i = data["i"]
    #     raan = data["RAAN"]
    #     omega = data["omega"]
    #     nu = data["nu"]

    #     self.visualisation_plane_4.setH(raan)
    #     self.visualisation_plane_4.setP(i-90)
        

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
        self.label_1 = self.add_text_label(text="label 1" , pos=(x_po , y_st))
        self.label_2 = self.add_text_label(text="label 2" , pos=(x_po , y_st-y_sp))
        self.label_3 = self.add_text_label(text="label 3" , pos=(x_po , y_st-y_sp*2))
        self.label_4 = self.add_text_label(text="label 4" , pos=(x_po , y_st-y_sp*3))
        self.label_5 = self.add_text_label(text="label 5" , pos=(x_po , y_st-y_sp*4))
        self.label_6 = self.add_text_label(text="label 6" , pos=(x_po , y_st-y_sp*5))
        self.label_7 = self.add_text_label(text="label 7" , pos=(x_po , y_st-y_sp*6))

        # Top right pause
        self.pause_label = self.add_text_label(text="II" , pos=(1.2 , y_st))
        self.pause_label.hide()
        self.time_label = self.add_text_label(text="_" , pos=(1 , y_st-y_sp))
        

        

        y_st = -0.9
        # Bottom left constants
        self.label_G = self.add_text_label(text=f"G = {G}" , pos=(x_po , y_st))
        self.label_M = self.add_text_label(text=f"M = {M}" , pos=(x_po , y_st+y_sp))
        self.label_J2 = self.add_text_label(text=f"J2 = {J2}" , pos=(x_po , y_st+y_sp*2))
        
        # Bottom right
        self.frames_left_label = self.add_text_label(text="_" , pos=(-x_po , y_st) , alignment_mode=TextNode.ARight)


    def update_hud(self):
        i = self.otv.elements.inclination
        raan = self.otv.elements.raan
        e = self.otv.elements.eccentricity
        arg_perigee = self.otv.elements.arg_perigee
        mean_anomaly = self.otv.elements.mean_anomaly
        a = self.otv.elements.a

        i = round(i , 3)
        raan = round(raan , 3)
        e = round(e , 3)
        arg_perigee = int(arg_perigee)
        mean_anomaly = int(mean_anomaly)
        a = round(a , 3)
        

        self.label_1.setText(f"i = {i}")
        self.label_2.setText(f"raan = {raan}")
        self.label_3.setText(f"e = {e}")
        self.label_4.setText(f"arg_perigee = {arg_perigee}")
        self.label_5.setText(f"mean_anomaly = {mean_anomaly}")
        self.label_6.setText(f"a = {a}")
        self.label_7.setText(f"")

        
        self.time_label.setText(f"t = {round(self.initial_delay,1)}")
        pc_done = int(100*len(self.log_a)/N_log)
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

    


if __name__ == "__main__":
    app = MyApp()
    app.run()