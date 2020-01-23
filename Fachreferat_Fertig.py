from __future__ import division
from visual import *
import math
from visual.graph import *

scene.width = 1200
scene.height = 450
scene.range = 14
scene.title = 'Rotating Box Object in Zero G'
scene.center = (5,0,0)
scene.background = color.white

box = box(pos=(-3,0,0), size=(2.0,3.5,1.0), omega=vector(4,0.1,0.1), mass=5.0, color=color.red)

Ixx = (box.mass/12)*(box.size.y*box.size.y+box.size.z*box.size.z)
Iyy = (box.mass/12)*(box.size.x*box.size.x+box.size.z*box.size.z)
Izz = (box.mass/12)*(box.size.x*box.size.x+box.size.y*box.size.y)

basis_vectors = []
basis_vectors.append(vector(1.,0,0))
basis_vectors.append(vector(0,1.,0))
basis_vectors.append(vector(0,0,1.))


xaxis = arrow(pos=(-3,0,0), axis=(2,0,0), shaftwidth=.1, color=color.yellow, opacity=.5)
yaxis = arrow(pos=(-3,0,0), axis=(0,2.75,0), shaftwidth=.1, color=color.green, opacity=.5)
zaxis = arrow(pos=(-3,0,0), axis=(0,0,1.5), shaftwidth=.1, color=color.blue, opacity=.5)

axis_List = []    
axis_List.append(xaxis)
axis_List.append(yaxis)
axis_List.append(zaxis)

xvec = arrow(pos=(10,0,0), length=(Ixx*box.omega.x)/5, axis=basis_vectors[0], shaftwidth=.1, color=color.yellow, opacity=.7)  
yvec = arrow(pos=(10+xvec.length,0,0), length=(2*Ixx*box.omega.y)/5, axis=basis_vectors[1], shaftwidth=.1, color=color.green, opacity=.7)
zvec = arrow(pos=(10+xvec.length,yvec.length,0), length=(Ixx*box.omega.z)/5, axis=basis_vectors[2], shaftwidth=.1, color=color.blue, opacity=.7)

totvec = arrow(pos=(10,0,0), length=sqrt(xvec.length*xvec.length+yvec.length*yvec.length+zvec.length*zvec.length), axis=norm((xvec.length,yvec.length,zvec.length)), shaftwidth=.2, color=color.orange)

vector_List = [] 
vector_List.append(xvec)
vector_List.append(yvec)
vector_List.append(zvec)
vector_List.append(totvec)


graph1 = gdisplay(title = 'Angular Momenta and Energy vs. Time',x=200,y=450,height = 325)

p_x = gcurve(gdisplay = graph1,color=color.yellow)
p_y = gcurve(gdisplay = graph1,color=color.green)
p_z = gcurve(gdisplay = graph1,color=color.blue)
p_tot = gcurve(gdisplay = graph1,color=color.orange)
E_tot = gcurve(gdisplay = graph1,color=color.white)

plot_List = [] 
plot_List.append(p_x)
plot_List.append(p_y)
plot_List.append(p_z)
plot_List.append(p_tot)
plot_List.append(E_tot)


def find_omega_dot(w):
    wx_dot = (Iyy-Izz)*w.y*w.z/Ixx
    wy_dot = (Izz-Ixx)*w.z*w.x/Iyy
    wz_dot = (Ixx-Iyy)*w.x*w.y/Izz
    w_dot = vector(wx_dot, wy_dot, wz_dot)
    return(w_dot)

def update_omega(w,wdot):
    w.x = w.x + dt*wdot.x
    w.y = w.y + dt*wdot.y
    w.z = w.z + dt*wdot.z

def rotate_box(obj,axes) :
    for i in range(len(axes)):
        obj.rotate(angle=obj.omega[i]*dt, axis=norm(axes[i].axis), origin=obj.pos)

def update_axis(List,w): #rotating the vectors and adjusting length to be proportional to angular momentum
    old = List
    for i in range(len(List)):
        k = 0
        while k<3:
            List[i].rotate(angle=w[k]*dt, axis=norm(old[k].axis), origin=List[i].pos)
            k += 1


def update_vectors(List,w):
    k = vector(0,0,0)
    for i in range(len(List)-1):
        k[i] = ((Ixx*w[i])/5)/abs((Ixx*w[i])/5)
        List[i].axis = k[i]*basis_vectors[i]
        List[i].length = abs((Ixx*w[i])/5)
        if i>0:
            List[i].pos.x = 10+(List[0].length)*k[0]
            if i>1:
                List[i].pos.y = (List[1].length)*k[1]
    List[3].axis = norm((k[0]*List[0].length, k[1]*List[1].length, k[2]*List[2].length))
    List[3].length = sqrt(List[0].length*List[0].length+List[1].length*List[1].length+List[2].length*List[2].length)
       

def update_plot(pList, w, time): 
    pList[0].plot(pos=(time, w.x*Ixx))
    pList[1].plot(pos=(time, w.y*Iyy))
    pList[2].plot(pos=(time, w.z*Izz))
    pList[3].plot(pos=(time, sqrt(w.x*w.x*Ixx*Ixx+w.y*w.y*Iyy*Iyy+w.z*w.z*Izz*Izz)))
    pList[4].plot(pos=(time, .5*(w.x*w.x*Ixx+w.y*w.y*Iyy+w.z*w.z*Izz)))

       
t = 0.0000
dt = .001

while t<60:
    rate(1000)
    omega_dot = find_omega_dot(box.omega)
    update_omega(box.omega, omega_dot)
    rotate_box(box,axis_List)
    update_axis(axis_List, box.omega)
    update_vectors(vector_List, box.omega)
    update_plot(plot_List, box.omega, t)
    t += dt
