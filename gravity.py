from math import floor, sqrt
import tkinter
from typing import Tuple
import customtkinter as ctk
import numpy as np
import time

class Planet:
    def __init__(self, mainwindow, xpos, ypos, radius, canvas, id):
        
        self.id = id
        self.mainwindow = mainwindow
        self.canvas = canvas
        self.x0 = xpos
        self.y0 = ypos
        self.radius = radius
        self.mass = radius**3
        self.r0 = np.array((xpos,ypos))
        self.velocity = np.array((0,0))
        self.thisplanet = self.canvas.create_aa_circle(self.x0, self.y0, self.radius)
        self.pinned = False

    def calculatespacegravity(self):
        if self.pinned:
            return
        force = np.array((0,0))

        for id, planet in enumerate(mainwindow.planets):
            if id != self.id:
                
                distance = np.linalg.norm(self.r0-planet.r0)
                direction = (planet.r0-self.r0)/distance
               
                force = force + direction * (self.mass)*(planet.radius**2)/(distance**2)
        self.velocity = self.velocity + 15*force/(self.mass)

    def calculatecursorgravity(self):
        if mainwindow.cursorx == -1:
            return
        cursor = np.array((mainwindow.cursorx, mainwindow.cursory))
        distance = max(np.linalg.norm(self.r0-cursor),self.radius)
        direction = (cursor-self.r0)/distance
        
        if mainwindow.walls_var.get() == "Wrap":
            for i in range(-1,2):
                for j in range(-1,2):
                    mirror = np.array((i*mainwindow.canvaswidth,j*mainwindow.canvasheight))
                    contenderdistance = max(np.linalg.norm(self.r0-(cursor+mirror)),self.radius)
                    if contenderdistance < distance:
                        distance = contenderdistance
                        direction = (cursor+mirror-self.r0)/distance

        

        force = direction * (self.mass) / (distance**1) / 10
        self.velocity = self.velocity + 3000*force / self.mass

    def calculateearthgravity(self):
        if self.pinned or self.r0[1]+self.radius>= mainwindow.canvasheight:
            return
        self.velocity[1] = self.velocity[1] + 3
        
    def calculatewalls(self):
        if (self.r0[0]-self.radius)<=0:
            self.velocity[0] = abs(self.velocity[0])*mainwindow.COR
        if (self.r0[0]+self.radius)>= mainwindow.canvaswidth:
            self.velocity[0] = -abs(self.velocity[0])*mainwindow.COR
        if (self.r0[1]-self.radius)<=0:
            self.velocity[1] = abs(self.velocity[1])*mainwindow.COR
        if (self.r0[1]+self.radius)>= mainwindow.canvasheight:
            self.velocity[1] = -abs(self.velocity[1])*mainwindow.COR
            if abs(self.velocity[1]) <1:
                self.velocity[0] = 0
                self.velocity[0] = self.velocity[0]*0.95
                
            
        if self.r0[0]-self.radius<0:
            self.r0[0] = self.radius
        elif self.r0[0]+self.radius>mainwindow.canvaswidth:
            self.r0[0] = mainwindow.canvaswidth - self.radius
        if self.r0[1]-self.radius<0:
            self.r0[1] = self.radius
        elif self.r0[1]+self.radius>mainwindow.canvasheight:
            self.r0[1] = mainwindow.canvasheight - self.radius
            
    def calculatewrap(self):
        if self.r0[0]<0:
            self.r0[0] = mainwindow.canvaswidth
        elif self.r0[0]>mainwindow.canvaswidth:
            self.r0[0] = 0
        if self.r0[1]<0:
            self.r0[1] = mainwindow.canvasheight
        elif self.r0[1]>mainwindow.canvasheight:
            self.r0[1] = 0
    
    def explode(self,epicenter):
        vec = self.r0-epicenter
        distance = np.linalg.norm(vec)
        direction = vec/distance
        self.velocity = self.velocity + direction/sqrt(distance) * 1000

    def moveplanet(self):
        if self.pinned:
            return
        self.r0 = self.r0 + 0.1*self.velocity
        wallversion = mainwindow.walls_var.get()
        if wallversion == "Walls": 
            self.calculatewalls()
        if wallversion == "Wrap":
            self.calculatewrap()
            
        self.canvas.coords(self.thisplanet, self.r0[0], self.r0[1])
        
    def radiuschange(self, radius):
        self.radius = radius
        self.canvas.coords(self.thisplanet, self.x0, self.y0, self.radius)

    def updatearrow(self,xpos,ypos):
        self.r0 = np.array((self.x0, self.y0))
        r1 = np.array((xpos,ypos))
        n = self.radius / 3
        self.rdir = self.r0-r1
        rdirlen = np.linalg.norm(self.rdir)
        rtip = self.r0 + self.rdir
        rperpendicular = np.array((-self.rdir[1],self.rdir[0]))/rdirlen
        
        p1 = self.r0+n*rperpendicular
        p2 = p1+self.rdir
        p3 = p2+n*rperpendicular
        p4 = p3 - 2*n*rperpendicular+4*n*self.rdir/rdirlen
        p5 = p3-4*n*rperpendicular
        p6 = p5+n*rperpendicular
        p7 = p6-self.rdir

        self.canvas.coords(self.thisarrow,
                              p1[0],p1[1],
                              p2[0],p2[1],
                              p3[0],p3[1],
                              p4[0],p4[1],
                              p5[0],p5[1],
                              p6[0],p6[1],
                              p7[0],p7[1])
    
        self.velocity = self.rdir/5

    def removearrow(self):
        self.canvas.delete(self.thisarrow)
        self.mass = self.radius**3

    def makearrow(self):
        self.thisarrow = self.canvas.create_polygon(self.x0,self.y0,self.x0,self.y0,self.x0,self.y0,self.x0,self.y0,self.x0,self.y0,self.x0,self.y0,self.x0,self.y0, fill = 'white')

    def pin(self):
        if not self.pinned:
            self.velocity = np.array((0,0))
            self.pinned = True
        else:
            self.pinned = False

class MainWindow(ctk.CTk):

    def drawcircle(self, event, x):
        radius = floor(x**2+x)+7
        
        def stoploop(loop):
            
            self.after_cancel(loop)
            self.maincanvas.unbind("<ButtonRelease-1")
        if x == 0:
            self.maincanvas.bind("<ButtonRelease-1>", lambda event: [stoploop(self.after_id), self.makearrow(-1,-1)])
            planet = Planet(self, event.x, event.y, 1, self.maincanvas, self.planetnum)
            self.planetnum += 1
            self.planets.append(planet)
        x += 0.3
        
        self.planets[-1].radiuschange(radius)
        if radius != -1:
            self.after_id = self.after(25, self.drawcircle, event, x)

    def makearrow(self, xpos, ypos):
        self.planets[-1].makearrow()
        self.maincanvas.unbind("<ButtonRelease-1>")
        self.maincanvas.bind("<Button-1>", lambda event: [self.maincanvas.unbind("<Motion>"), 
                                                          self.maincanvas.bind("<Button-1>", lambda event: self.drawcircle(event, 0)),
                                                          self.planets[-1].removearrow()])
        self.maincanvas.bind("<Motion>", lambda event: self.planets[-1].updatearrow(event.x,event.y))

    def explosion(self,event):
        for planet in self.planets:
            planet.explode(np.array((event.x,event.y)))

    def togglegame(self):
        def stoploop(loop):
            try:
                self.after_cancel(loop)
            except ValueError:
           
                return
        if self.toggle == 0:
            self.toggle = 1
            self.startbtn.configure(text = "Stop")
            self.version.configure(state="disabled")
            self.CORentry.configure(state = "disabled")
            self.gridentry.configure(state = "disabled")
            self.maincanvas.bind("<Button-1>", self.explosion)
            if self.version_var.get() == "Cursor Gravity":
                self.maincanvas.bind("<Motion>", self.getcursorpos)
                self.cursorx = -1
                self.cursory = -1
            
            self.rungame()
        elif self.toggle == 1:
            self.toggle = 0
            self.startbtn.configure(text = "Run")
            self.version.configure(state="normal")
            self.CORentry.configure(state = "normal")
            self.gridentry.configure(state = "normal")
            self.maincanvas.unbind("<Motion>")
            self.maincanvas.bind("<Button-1>", lambda event: self.drawcircle(event, 0))
            stoploop(self.afterloop)

    def createchuncks(self):
        self.chunks = []
        packingfactor = 1
        
        numofchunks = max((self.planetnum)//packingfactor,1)
        (width, height) = self.getrect(numofchunks)
        for i in range(int(width)):
            self.chunks.append([])
            for j in range(int(height)):
                self.chunks[i].append([])
        self.chunkingwidth = self.canvaswidth / width
        self.chunkingheight = self.canvasheight / height
        for planet in self.planets:
            planet.entergrid()
                
    def getcursorpos(self, event):
        self.cursorx = event.x
        self.cursory = event.y

    def rungame(self):
        self.simulate()
        self.calculateforces()
        for planet in self.planets:
            planet.moveplanet()
        self.afterloop = self.after(10, self.rungame)

    def simulate(self):
        for n, planet1 in enumerate(self.planets):
            for planet2 in self.planets[n+1::]:
                distance = np.linalg.norm(planet1.r0-planet2.r0)
                if distance <= planet1.radius+planet2.radius:
                    self.collision(planet1, planet2)
        
    def calculateforces(self):
        for planet in self.planets:
            if self.version_var.get() == "Space Gravity":
                planet.calculatespacegravity()
                
            elif self.version_var.get() == "Earth Gravity":
                planet.calculateearthgravity()
                
            elif self.version_var.get() == "Cursor Gravity":
                planet.calculatecursorgravity()      

    def collision(self,a,b):
        vec_ab = a.r0-b.r0
        overlap = a.radius + b.radius - np.linalg.norm(vec_ab)
        if np.linalg.norm(vec_ab) == 0:
            return
        loi = vec_ab/np.linalg.norm(vec_ab)
        normal = np.array((-loi[1],loi[0]))
        v_ai = np.dot(a.velocity,loi)*loi
        v_bi = np.dot(b.velocity,loi)*loi
        v_an = np.dot(a.velocity,normal)*normal
        v_bn = np.dot(b.velocity,normal)*normal
        if not a.pinned:
            v_af = (v_ai*a.mass+v_bi*b.mass+b.mass*self.COR*(v_bi-v_ai))/(a.mass+b.mass)
            a.velocity = v_af+v_an
            a.r0 = a.r0 + overlap/2 * loi
        if not b.pinned:
            v_bf = (v_bi*b.mass+v_ai*a.mass+a.mass*self.COR*(v_ai-v_bi))/(a.mass+b.mass)
            b.velocity = v_bf+v_bn
            b.r0 = b.r0 - overlap/2 * loi
        
    def clear(self):
        self.planetnum = 0
        self.planets = []
        self.maincanvas.delete("all")

    def pinplanet(self, event): 
        for planet in self.planets:
            if np.linalg.norm(planet.r0-np.array((event.x,event.y)))<planet.radius:
                planet.pin()

    def CORentered(self, event):
        try:
            self.COR = float(self.CORentry.get())
        except ValueError:
            print("Please enter a number")
            self.CORentry.delete(0, ctk.END)
        else:
            self.focus()

    def getrect(self, n):
        n = int(n)
        factors = [i for i in range(1,n+1) if n%i==0]
        if len(factors)%2 == 1:
            height = sqrt(n)
            width = sqrt(n)
        else:
            differences = [factors[-i-1]-factors[i] for i in range(len(factors)//2)]
            height = int(factors[differences.index(min(differences))])
            width = int(height + min(differences))
        return width, height

    def gridentered(self, event):
        try:
            n = int(self.gridentry.get())
            (width, height) = self.getrect(n)
        except ValueError:
            print("Please enter an integer")
            self.gridentry.delete(0, ctk.END)
            return
        self.clear()
        self.focus()
        r = 0.5
        xspacing = self.canvaswidth*r/(width-1)
        if height != 1:
            yspacing = self.canvasheight*r/(height-1)
            radius = int(min(xspacing,yspacing)//6)
        else:
            radius = int(xspacing//6)
        for i in range(n):
            x = self.canvaswidth*(1/2-r/2) + i%width * xspacing
            if height != 1:
                y = self.canvasheight*(1/2-r/2) + i//width * yspacing
            elif height == 1:
                y = self.canvasheight/2
            planet = Planet(self, x, y, radius, self.maincanvas, self.planetnum)
            self.planets.append(planet)
            self.planetnum += 1
            
    def configure(self,event):
        
        self.canvaswidth = event.width
        self.canvasheight = event.height
     

    def initialize(self):
   
        self.btnframe = ctk.CTkFrame(self)
        self.columnconfigure(0, weight = 1)
        self.columnconfigure(1, weight = 1)
        self.columnconfigure(2, weight = 1)
        self.rowconfigure(0, weight = 1)
        self.rowconfigure(1, weight = 1)
        self.startbtn=ctk.CTkButton(self.btnframe, text="Run", command = self.togglegame)
        self.startbtn.grid(column = 0, row = 0, padx = 20, pady = 20)
        self.clearbtn=ctk.CTkButton(self.btnframe, text = "Clear", command = self.clear)
        self.clearbtn.grid(column = 1, row = 0, padx = 20, pady = 20)
        self.version_var = ctk.StringVar(value = "Space Gravity")
        self.version = ctk.CTkOptionMenu(self.btnframe, values = ["Space Gravity", "Earth Gravity", "Cursor Gravity", "No Gravity"], variable = self.version_var)
        self.version.grid(column = 1, row = 1, padx = 20, pady = 20)
        self.walls_var = ctk.StringVar(value = "Walls")
        self.wallcheck = ctk.CTkOptionMenu(self.btnframe, values = ["Walls","Wrap","None"], variable = self.walls_var)
        self.wallcheck.grid(column = 0, row = 1, padx = 20, pady = 20)
        self.CORentry = ctk.CTkEntry(self.btnframe, placeholder_text="Coefficient of Restitution", width=160)
        self.CORentry.grid(column = 2, row = 0, padx = 20, pady = 20)
        self.CORentry.bind("<Return>", self.CORentered)
        self.gridentry = ctk.CTkEntry(self.btnframe, placeholder_text="Create Grid of Balls", width=160)
        self.gridentry.grid(column = 2, row = 1, padx = 20, pady = 20)
        self.gridentry.bind("<Return>", self.gridentered)
        self.btnframe.pack(pady=10)
        self.maincanvas=ctk.CTkCanvas(self, width = self.canvaswidth, height = self.canvasheight, bg='#181818')
        self.maincanvas.pack(fill = ctk.BOTH, expand = True)
        self.maincanvas.bind("<Button-1>", lambda event: self.drawcircle(event, 0))
        self.maincanvas.bind("<Button-3>", self.pinplanet)
        self.maincanvas.bind("<Configure>",self.configure)
        self.planets = []
        self.planetnum = 0
        self.toggle = 0
        self.afterloop = None
        self.COR = 0.8
        self.collisionshappening = []

    def __init__(self, fg_color: str | Tuple[str, str] | None = None, **kwargs):
        super().__init__(fg_color, **kwargs)
        screen_width = 1000
        screen_height = 800
        self.canvaswidth = screen_width-200
        self.canvasheight = screen_height-200
        self.geometry(f"{screen_width}x{screen_height}+0+0")
        self.title("Gravity Simulation")
        self.firstrun = True
        self.initialize()
        
    
mainwindow = MainWindow()
mainwindow.mainloop()