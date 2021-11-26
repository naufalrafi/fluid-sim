import Python3.myClasses as mc

#Setting up boundary conditions for horizontal velocity
def setUBoundary(space: mc.Space, left: mc.Boundary, right: mc.Boundary, top: mc.Boundary, bottom: mc.Boundary):
    if left.type == 'D':
        space.u[:,0] = left.value
    elif left.type == 'N':
        space.u[:,0] = (-left.value * space.dx) + space.u[:,1]

    if right.type == 'D':
        space.u[:,-1] = right.value
    elif right.type == 'N':
        space.u[:,-1] = (right.value * space.dx) + space.u[:,-2]

    if top.type == 'D':
        space.u[-1,:] = (2 * top.value) - space.u[-2,:]
    elif top.type == 'N':
        space.u[-1,:] = (-top.value * space.dy) + space.u[-2,:]
     
    if bottom.type == 'D':
        space.u[0,:] = (2 * bottom.value) - space.u[1,:]
    elif bottom.type == 'N':
        space.u[0,:] = (bottom.value*space.dy) + space.u[1,:] 


#Setting up boundary conditions for verticaal velocity
def setVBoundary(space: mc.Space, left: mc.Boundary, right: mc.Boundary, top: mc.Boundary, bottom: mc.Boundary):
    if left.type == 'D':
        space.v[:,0] = (2 * left.value) - space.v[:,1]
    elif left.type == 'N':
        space.v[:,0] = (-left.value * space.dx) + space.v[:,1]
    
    if right.type == 'D':
        space.v[:,-1] = (2 * right.value) - space.v[:,-2]
    elif right.type =='N':
        space.v[:,-1] = (right.value*space.dx) + space.v[:,-2]
        
    if top.type == 'D':
        space.v[-1,:] = top.value
    elif top.type == 'N':
        space.v[-1,:] = (-top.value * space.dy) + space.v[-2,:]
     
    if bottom.type == 'D':
        space.v[0,:] = bottom.value
    elif bottom.type=="N":
        space.v[0,:] = (bottom.value * space.dy) +space.v[1,:]

#Seting up boundary conditions for pressure
def setPBoundary(space: mc.Space, left: mc.Boundary, right: mc.Boundary, top: mc.Boundary, bottom: mc.Boundary):
    if left.type == 'D':
        space.p[:,0] = left.value
    elif left.type == 'N':
        space.p[:,0] = (-left.value * space.dx) + space.p[:,1]
    
    if right.type == 'D':
        space.p[1,-1] = right.value
    elif right.type == 'N':
        space.p[:,-1] = (right.value * space.dx) + space.p[:,-2]
        
    if top.type == 'D':
        space.p[-1,:] = top.value
    elif top.type == 'N':
        space.p[-1,:] = (-top.value * space.dy) + space.p[-2,:]
     
    if bottom.type == 'D':
        space.p[0,:] = bottom.value
    elif bottom.type == 'N':
        space.p[0,:] = (bottom.value * space.dy) + space.p[1,:]