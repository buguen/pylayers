#!/usr/bin/python
#-*- coding:Utf-8 -*-
from serial import Serial
import pdb
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class Profile(object):
    Accmax = 200
    Vmax = 10
    Dmax = 1500

    def __init__(self,**kwargs):
        """

        Parameters
        ----------

        num     : profile number
        aa      : aceleration (rps2)
        ad      : deceleration (rps2)
        dstep   : distance (step)
        vmax    : vmax (rps)
        vs      : vstart (rps)
        spr     : steps per round

        """

        defaults = {'num': 1,
                    'aa': 200,
                    'ad': 200,
                    'dstep': 12800,
                    'vmax': 15,
                    'vs': 0,
                    'spr': 4000,
                    'N':100}

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        self.num = kwargs['num']
        self.aa = kwargs['aa']
        self.ad = kwargs['ad']
        self.dstep = kwargs['dstep']
        self.vmax = kwargs['vmax']
        self.vs = kwargs['vs']
        self.spr = kwargs['spr']
        self.N = kwargs['N']

        #
        # spr   : steps per revolution
        # dstep : number of steps
        # drev  : number of revolutions
        #

        self.drev = self.dstep/(1.0*self.spr)
        self.T = (self.drev+self.vmax**2/self.aa)/(1.0*self.vmax)

        #assert(0<accmax)
        #assert(0<v<vmax)
        #assert(0<d<dmax)

        self.cmd = 'PROFILE'+str(self.num)+\
                '('+str(self.aa)+','+str(self.ad)+\
                ','+str(self.dstep)+','+str(self.vmax)+\
                ','+str(self.vs)+')'

        self.t = np.linspace(0,self.T,self.N)
        self.v = self.vmax*np.ones(self.N)
        self.v1 = self.t*self.aa
        self.v2 = self.t*self.ad

        #
        # t1 : end of acceleration phase
        # t2 : begining of deceleration phase
        #
        t1 = self.vmax/(1.0*self.aa)
        t2 = self.T-self.vmax/(1.0*self.ad)

        u1 = np.where(self.t<t1)[0]
        u2 = np.where(self.t>=t2)[0]

        self.v[u1] = self.t[u1]*self.aa
        self.v[u2] = -self.t[u2]*self.ad+(self.vmax+self.ad*t2)
        #
        # dr : distance in number of revolutions
        # ds : distance in number of steps
        #
        self.dr = np.cumsum(self.v)*(self.t[1]-self.t[0])
        self.ds = self.dr*self.spr

    def __repr__(self):
        st = ''
        st = st + 'num : '+ str(self.num) + '\n'
        st = st + 'acceleration  : '+ str(self.aa) + '\n'
        st = st + 'deceleration  : '+ str(self.ad) + '\n'
        st = st + 'dstep : '+ str(self.dstep) + '/'+ str(self.ds[-1]) + '\n'
        st = st + 'maximum velocity : '+ str(self.vmax) + '\n'
        st = st + 'starting velocity : '+ str(self.vs) + '\n'
        st = st + 'spr : '+ str(self.spr) + '\n'
        st = st + 'N  : '+ str(self.N) + '\n'
        return(st)

    def duration(self):
        """
        """
        self.duration

    def show(self):
        """Enables view profile

        Examples
        --------

        P1=Profile(1,200,200,15000,15,0)
        P1.show()

        """
        plt.subplot(211)
        plt.plot(self.t,self.v)
        plt.xlabel('time (s)')
        plt.ylabel('Velocity (rev/s)')
        plt.title('Evolution of velocity w.r.t time')
        plt.subplot(212)
        plt.plot(self.t,self.ds)
        plt.xlabel('time (s)')
        plt.ylabel('distance (step)')
        plt.title('Evolution of distance w.r.t time')
        plt.legend()
        plt.show()

class Axes(object):
    svar  = {'BU':'Buffer Usage',
            'CQ':'Command queuing',
            'DF':'Drive Fault status',
            'EI':'Encoder Input',
            'EO':'Encoder signal output',
            'EP':'Encoder Position',
            'ER':'Feedback encoder resolution',
            'EX':'Coms response Style & echo control',
            'IN':'Inputs',
            'IP':'In position flag',
            'IT':'IN position Time',
            'MC':'Motor current',
            'MR':'Motor resolution',
            'MS':'Motor Standby current',
            'MV':'Moving',
            'PA':'Position Absolute',
            'PE':'Position Error',
            'PI':'Position Incremental',
            'RB':'Ready/Busy flag',
            'RM':'Registration Move',
            'RV':'Revision software',
            'SN':'Serial Number',
            'ST':'Status of indexing',
            'UF':'User Program Fault status'
           }


    dstatus = {}
    dstatus[1]='command processing paused'
    dstatus[2]='looping'
    dstatus[3]='wait for trigger'
    dstatus[4]='running program'
    dstatus[5]='going home'
    dstatus[6]='waiting for delay timeout'
    dstatus[7]='registration in progress'
    dstatus[9]='motor energised'
    dstatus[11]='event trigger active until trigger inputs are reset'
    dstatus[12]='input in LSEL not matching label'
    dstatus[13]='-ve limit seen during last move'
    dstatus[14]='+ve limit seen during last move'
    dstatus[19]='moving'
    dstatus[20]='stationnary'
    dstatus[21]='no registration signal seen in registration window'
    dstatus[22]='cannot stop within the defined regisration distance'

    dusrflt = {}
    dusrflt[1]='value is out of range'
    dusrflt[2]='incorrect command syntax'
    dusrflt[3]='last label already in use'
    dusrflt[4]='label of this name not defined'
    dusrflt[5]='missing Z pulse when homing'
    dusrflt[6]='homing failed no signal detected'
    dusrflt[7]='home signal too narrow'
    dusrflt[8]='drive de-energised'
    dusrflt[9]='cannot related end statement to a label'
    dusrflt[10]='program memory buffer full'
    dusrflt[11]='no more motion profiles avaible'
    dusrflt[12]='no more sequence labels avaible'
    dusrflt[13]='end of travel limit hit'
    dusrflt[14]='still moving'
    dusrflt[15]='deceleration error'
    dusrflt[16]='transmit buffer overflow'
    dusrflt[17]='user program nesting overflow'
    dusrflt[18]='cannot use an undefined profile'
    dusrflt[19]='drive not ready'
    dusrflt[22]='save error'
    dusrflt[23]='command not supported by this product'

    ddrvflt={}
    ddrvflt[1]='Composite Fault'
    ddrvflt[2]='Output stage over curent'
    ddrvflt[3]='Supply rail failure'
    ddrvflt[4]='Ambient over temperature'
    ddrvflt[5]='Drive over temperature'
    ddrvflt[6]='Configuration error'
    ddrvflt[7]='Motor high voltage rail failure'
    ddrvflt[8]='Output fault'

    def __init__(self,_id,name,ser,scale=12800,typ='t'):
        """
        _id  : axes id
        name : axes name
        ser  : serial port
        scale : nstep/cm if typ='t'
        typ : 't'|'r'

        """
        self.status = [0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0]
        self.usrflt = [0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0]
        self.drvflt = [0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0]
        self._id = _id
        self.name = name
        # serial port
        self.ser = ser
        # list of profiles
        self.scale = scale
        self.typ = typ
        self.lprofile=[]


    def __repr__(self):
        st = 'st'
        st = st+str(self._id)
        return(st)

    def show(self):
        """
        """
        pass

    #def com2(self,com='R(ST)'):
        #self.ser.write(str(self._id)+com+'\r\n')
        #st = self.ser.readlines()
        #return(st)

    def getvar(self,lvar=[]):
        """
        Parameters
        ----------

        lvar : list of variables

        """
        if lvar == []:
            lvar = Axes.svar.keys()
        for var in lvar:
            st = self.com('R('+var+')')
            print Axes.svar[var],st[1]


    def com(self,command='R(SN)',verbose=False):
        """ send command to serial port

        Parameters
        ----------

        prefix : str command prefix
        arg :  command argument
        verbose :

        """
        cst = str(self._id) + command + '\r\n'
        if verbose:
            print cst
        self.ser.write(cst)
        st = self.ser.readlines()
        return(st)



    def home(self,cmd='get',**kwargs):
        """ enables back home

        Parameters
        ----------

        cmd  : 'get','set','go'

        """
        defaults = {'mode':0,
                    'vel':10,
                    'acc':10,
                    'edg':'+',
                    'typ':0,
                    'armed':1
        }

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        if cmd=='get':
            st = self.com('HOME')
            ans = st[1].split(' ')

            if '1' in ans[0]:
                print "armed, "
            else:
                print "not armed, "
            if '-' in ans[1]:
                print "reference edge is negative, "
            else:
                print "reference edge is positive, "
            if '1' in ans[2]:
                print "home switch normally closed 1, "
            else:
                print "home switch normally open 0 (default),  "

            if '+' in ans[3]:
                print 'velocity : +',eval(ans[3].split('V+')[1]), "rps"
                #print ", "
            else:
                print 'velocity : -',eval(ans[3].split('V-')[1]), "rps"

            print 'acceleraton : ',eval(ans[4].split('A')[1]), "rps²"

            if '0' in ans[5]:
                print 'Mode 0: Motor in the active window of the switch(default)'
            if '1' in ans[5]:
                print 'Mode 1: Motor in the position to the edge + or -'
            if '2' in ans[5]:
                print 'Mode 2: Improve homing repeatability'

        if cmd=='set':

            if kwargs['vel']>0:
                vel = '+'+str(kwargs['vel'])
            else:
                vel = '-'+str(kwargs['vel'])


            cstr = 'HOME'+str(kwargs['armed'])+\
                          '('+kwargs['edg']+','+\
                          str(kwargs['typ'])+','+\
                          vel+','+\
                          str(kwargs['acc'])+','+\
                          str(kwargs['mode'])+')'

            self.com(cstr)

        if cmd=='go':
            cstr = 'HOME1'
            self.com(cstr)
            cstr = 'ARM1'
            self.com(cstr)
            cstr = 'GH'
            self.com(cstr)

    def del_profile(self,index=1):
        """ delete profile

        Parameters
        ----------

        index : int

        """
        self.lprofile.pop(index-1)

    def add_profile(self,**kwargs):
        """  add a new profile to lprofile

        Parameters
        ----------

        num     : profile number
        aa      : aceleration (rps2)
        ad      : deceleration (rps2)
        dstep   : distance (step)
        v       : vmax (rps)
        vs      : vstart (rps)
        spr     : steps per round

        """
        defaults = {'num': 0,
                    'aa': 200,
                    'ad': 200,
                    'dstep': 12800,
                    'v': 15,
                    'vs': 0,
                    'spr': 4000,
                    'N':100}

        for k in defaults:
            if k not in kwargs:
                kwargs[k]=defaults[k]

        num = len(self.lprofile)
        kwargs['num']=num+1
        if len(self.lprofile)<8:
            prof = Profile(**kwargs)
            # update profile
            self.lprofile.append(prof)
            # send command 
            # self.com(prof.cmd)

    def set_profile(self,num):
        """
        """
        assert(num<=len(self.lprofile)),"profile number not defined" 
        self.com(str(self._id)+'USE'+str(num))

    def reset(self):
        """ reset axis
        """
        self.com('OFF')
        self.com('ON')



    def mvpro(self,id_pro):
        """ move according to specified profile

        Parameters
        ----------

        axis
        id_pro

        Returns
        -------

        str

        Examples
        --------

        sm.mvpro(1,1)  # axis 1 , profile 1
        """

        #com = 'USE('+str(id_pro)+')'
        #self.com(com)
        #self.com('G')
        
        com = 'USE('+str(id_pro)+')'
        self.com(com)
        self.com('G')

    def read(self):
        pass

    def reg(self,typ='ST'):
        """ read boolean quantities in registers  : ST,UF,DF
    
        Examples
        --------
          
        s.a[1].reg('ST') #scans over axis 1 by given status
        """

        buf = self.com('R','('+typ+')')
        buf = buf[1]
        buf = buf.replace('*','').replace('\r\n','').split('_')

        for k in range(8):
            for l in range(4):
                val = eval(buf[k][l])
                if typ=='ST':
                    self.status[k*4+l] = val
                    if val:
                        print Axes.dstatus[k*4+l+1]
                if typ=='UF':
                    self.usrflt[k*4+l] = val
                    if val:
                        print Axes.dusrflt[k*4+l+1]
                if typ=='DF':
                    self.drvflt[k*4+l] = val
                    if val:
                        print Axes.ddrvflt[k*4+l+1]

    def mv(self,var=0):
        """ move axes in translation or rotation
        
        Parameters
        ----------
        
        var : distance (cm) | degres (°)
       
        Examples
        --------

        s.a[1].mv(10)  # moves over 10cm on axis 1
        s.a[3].mv(45) # moves over 45° on axis 3
        """
        #assert(self.typ=='t'),'Axes is not a linear axes'
        #nstep = dcm*self.scale
        #com = self.com('D'+str(nstep))
        #if typ=='t':
            #nstep = dcm*self.scale
            #com = self.com('D'+str(nstep))
        
        nstep = int(var*self.scale)
        scom1 = 'D'+str(nstep)
        com = self.com(scom1)
        com = self.com('G')
        #com = self.com(scom1,verbose=True)
        #scom2 = 'G'es
        #com = self.com(scom2,verbose=True)
        #print "distance parcourue : ", var+str('cm')  
        #com = self.com(scom1,verbose=True)
        com = self.com(scom2,verbose=True)
        com = self.com(scom3,verbose=True)
    def homeor(self):
        """Back to material origin
        """
        scom0 = 'HOME'+ str(self._id)
        scom1 = 'ARM'+ str(self._id)
        com = self.com(scom0)
        com = self.com(scom1)
        com = self.com('GH')

    def homing(self,typ=''):
        """Set up PA
        """
        pass

        #dstatus[13]='-ve limit seen during last move'
        #dstatus[14]='+ve limit seen during last move'
        #while Axes.dstatus[13]!= 1:
        #lire le buffer
        #buf = self.com('R','('+typ+')')
        #buf = buf[1]
        #buf = buf.replace('*','').replace('\r\n','').split('_')
        #print buf
        #if Axes.dstatus[13]== 0:
            #if Axes.dstatus[14]== 1:
                #scom1 = 'D-793600'
                #scom2 = 'G'
                #com   = self.com('D-793600',verbose=True)
                #com   = self.com('GH',verbose=True)
                #com   = self.com('R(PA)',verbose=True)
            #if :
        #else:
            #com   = self.com('R(PA)',verbose=True)

    def close(self):
        self.ser.close()


    def fromfile(self,cmdfile,dirfile='./DriverFiles'):
        cmdfile = dirfile+'/'+cmdfile
        fd = open(cmdfile,'r')
        lis = fd.readlines()
        for li in lis:
            li=li.replace("\n","\r\n")
            print li
            self.ser.write(li)
        st = self.ser.read(100)
        #while self.ser.inWaiting()>0:
        #    st += self.ser.read(10)
        #self.ser.close()
        return(st)
    #def translation(self,axis,offset):

class Scanner(object):
    def __init__(self,port):
        self.ser = Serial(port = port, baudrate=9600, timeout = 1)
        self.a  = ['',Axes(1,'x',self.ser,scale=12800),
                      Axes(2,'y',self.ser,scale=22800),
                      Axes(3,'rot',self.ser,scale=2111.1111111111113,typ='r')] #self.a4  = Axes(4,'z',self.ser,typ='r')

    def home(self):
        pass



if __name__=="__main__":
    pass
    #s = Scanner('/dev/ttyUSB2')
    #s = Scanner('/dev/ttyUSB1')
    #sm.fromfile('prog1')
    #sm.fromfile('AY')
    #Sc[1].com('ON')
    #st = sm.com(1,'LIMITS',(0,1,1))
    #st = sm.com(1,'1D-4000')
    #st = sm.com(1,'G')
