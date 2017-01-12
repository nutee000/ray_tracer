import numpy as np
from ray import *


class Plane:
    def __init__(self,n,D,ref,mn,color,Ka=[0.0,0.0,0.0],Ks=[0.5,0.5,0.5],Kd=[0.7,0.7,0.7],Ksn=7):
        self.n=self.norm(np.array(n))
        self.D=D
        self.ref=ref #ref=0 reflect, ref=1 refract
        self.color=np.array(color)
        self.oid=None
        self.mn=mn
        self.Ka=np.array(Ka)
        self.Ks=np.array(Ks)
        self.Kd=np.array(Kd)
        self.Ksn=Ksn
    def cal_inters(self,ray):
        temp=np.sum(self.n*ray.Rt)
        if temp==0:
            return None
        else:
            t=-(self.D+np.sum(self.n*ray.R0))/np.sum(self.n*ray.Rt)
            if t<=0.001:
                return None    #[refle0,3d1,n2,rw3,rray4,obj5,t6,rfray7]
            else:
                inters3d=ray.R0+t*ray.Rt
                rray=Ray()
                rfray=Ray()
                rfray.Rt=ray.Rt-np.sum(ray.Rt*self.n)*self.n*2
                rfray.R0=inters3d                
                if self.ref==0:
                    rray=np.copy(rfray)               
                else:
                    cost=(1-(1-(np.sum(self.n*ray.Rt))**2)/(self.mn**2))**0.5
                    if np.sum(ray.Rt*self.n)<=0:
                        rray.Rt=ray.Rt/self.mn-(np.sum(ray.Rt*self.n)/self.mn+cost)*self.n
                    else:
                        rray.Rt=ray.Rt/self.mn-(np.sum(ray.Rt*self.n)/self.mn-cost)*self.n
                    rray.R0=inters3d
                    print "obj",rray.R0,rray.Rt
                return [self.ref,inters3d,self.n,0.5,rray,self.oid,t,rfray]
    def norm(self,n):
        nsum=np.sqrt(np.sum(n*n))
        return n/nsum


class Sphere:
    def __init__(self,P,r,ref,mn,color,Ka=[0.0,0.0,0.0],Ks=[0.5,0.5,0.5],Kd=[0.7,0.7,0.7],Ksn=7):
        self.P=np.array(P)
        self.r=r
        self.ref=ref #ref=0 reflect, ref=1 refract
        self.color=np.array(color)
        self.oid=None
        self.mn=mn
        self.Ka=np.array(Ka)
        self.Ks=np.array(Ks)
        self.Kd=np.array(Kd)
        self.Ksn=Ksn
    def cal_inters(self,ray):
        l=self.P-ray.R0
        tp=np.sum(l*ray.Rt)
        if tp<=0 and np.sqrt(np.sum(l*l))+0.001>=self.r:
            return None
        else:
            d=np.sqrt(np.sum(l*l)-tp*tp)
            if d+0.001>=self.r:
                return None
            else:
                t1=np.sqrt(self.r*self.r-d*d)
                t=None
                if np.sqrt(np.sum(l*l))>=self.r+0.001:
                    t=tp-t1
                else:
                    t=t1+tp                
                inters3d=ray.R0+t*ray.Rt
                n=self.norm(inters3d-self.P)
                rray=Ray()
                rfray=Ray()
                rfray.Rt=ray.Rt-np.sum(ray.Rt*n)*n*2
                rfray.R0=inters3d                
                if self.ref==0:
                    rray=np.copy(rfray)
                else:
                    cost=(1-(1-(np.sum(n*ray.Rt))**2)/(self.mn**2))**0.5
                    if np.sum(ray.Rt*n)<=0:
                        rray.Rt=ray.Rt/self.mn-(np.sum(ray.Rt*n)/self.mn+cost)*n
                    else:
                        rray.Rt=ray.Rt/self.mn-(np.sum(ray.Rt*n)/self.mn-cost)*n
                    rray.R0=inters3d
                return [self.ref,inters3d,n,0.5,rray,self.oid,t,rfray]                    
    def norm(self,n):
        nsum=np.sqrt(np.sum(n*n))
        return n/nsum

class Tri:
    def __init__(self,v0,v1,v2,ref,mn,color,Ka=[0.0,0.0,0.0],Ks=[0.5,0.5,0.5],Kd=[0.7,0.7,0.7],Ksn=7):
        self.v=[np.array(v0),np.array(v1),np.array(v2)]
        tn=np.cross(self.v[1]-self.v[0],self.v[2]-self.v[1])
        self.n=self.norm(tn)
        self.ref=ref #ref=0 reflect, ref=1 refract
        self.color=np.array(color)
        self.oid=None
        self.mn=mn
        self.Ka=np.array(Ka)
        self.Ks=np.array(Ks)
        self.Kd=np.array(Kd)
        self.Ksn=Ksn
    def cal_inters(self,ray):
        e1=self.v[0]-self.v[1]
        e2=self.v[0]-self.v[2]
        s= self.v[0]-ray.R0
        a0=np.linalg.det(np.vstack((ray.Rt,e1,e2)))
        if a0==0.0:
            return None
        else:
            t=np.linalg.det(np.vstack((s,e1,e2)))/a0
            b=np.linalg.det(np.vstack((ray.Rt,s,e2)))/a0
            c=np.linalg.det(np.vstack((ray.Rt,e1,s)))/a0
            if t>0.001 and b>0 and b<1 and c>0 and c<1 and b+c<1:#[refle0,3d1,n2,rw3,rray4,obj5,t6,rfray7]
                inters3d=ray.R0+t*ray.Rt
                rray=Ray()
                rfray=Ray()
                rfray.Rt=ray.Rt-np.sum(ray.Rt*self.n)*self.n*2
                rfray.R0=inters3d                
                if self.ref==0:
                    rray=np.copy(rfray)
                else:
                    cost=(1-(1-(np.sum(self.n*ray.Rt))**2)/(self.mn**2))**0.5
                    if np.sum(ray.Rt*self.n)<=0:
                        rray.Rt=ray.Rt/self.mn-(np.sum(ray.Rt*self.n)/self.mn+cost)*self.n
                    else:
                        rray.Rt=ray.Rt/self.mn-(np.sum(ray.Rt*self.n)/self.mn-cost)*self.n
                    rray.R0=inters3d
                return [self.ref,inters3d,self.n,0.5,rray,self.oid,t,rfray]
            else:
                return None
    def norm(self,n):
        nsum=np.sqrt(np.sum(n*n))
        return n/nsum



if __name__=="__main__":
    print "hello nut"