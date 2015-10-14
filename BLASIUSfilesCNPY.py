import sys
sys.path.insert(0,'/home/s1359318/rls_datastore/PhD/BLASIUS/blasius_v4.01/2-d-canopy/Output')

# Above is the file location of Blasiusfiles.py
# This needs adding to the header of every python script

from matplotlib import rcParams
    
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator 

from mpl_toolkits.axes_grid1.inset_locator import* 

import numpy as np

#--------------------------------------------------------------------#

def openBlasiusfile(fname):
#def openBlasiusfile(fname, endian):
    
#  *endian code needs testing out
    
#    if ('little' in endian ):
#        print 'endian is little'
#    elif ('big' in endian):
#        print 'endian is big'

#    print args

        import numpy as np
    
        global numitems,fid
    
        fid=open(fname)
        
# If there are more than 200 items in the header this will need increasing. Better doing this rather than .append on each new entry
# the header will be saved as an array (see just below) consisting of name, offset, ndim and dim. dim covers the last 3 colums rather than 1 column where each entry has 3 dimensions. If the later is done then each entry is identicle. 

        file = np.array(np.zeros((200,6), dtype=list))

        numitems=0
        offset=0
        ndim=0
        
        dim = np.zeros(3) #creates an array with one row and three columns, with zeros in each column      
        
        dim[0]=1
        dim[1]=1
        dim[2]=1
        
        inchar = fid.read(1) #assigns one (the first) letter (!) in the header to the variable inchar, if you were to run inchar=fid.read(10) it would assign the first 10 letters etc
       
        tmpstr = np.array(np.zeros(200), dtype=str)   #creates an empty array of zeroes, with string data type
        

        while inchar != ':':
            
            if inchar == '!':
                
                #rosacount += 1            
                #print rosacount
            
                inchar = fid.read(1)
                #print inchar
                while inchar != '!':
                    inchar = fid.read(1)
                    #print inchar
                inchar=fid.read(1)
                #print inchar
                #print 'done this' 
                    
                           
                
            elif inchar == '&':
                
                #rosacount += 1            
                #print rosacount
                
                ndim=0
                dim[0]=1
                dim[1]=1
                dim[2]=1
                ccount=0
                inchar=fid.read(1)
                #print inchar
                
                
            elif ((isinstance(inchar, str)) and (inchar >= '0') and (inchar<='9')):

                #rosacount += 1            
                #print rosacount                
                
                if ndim != 0:
                    
                    print 'Error reading header'
                    fid.close()    
                    
                while ((inchar==',') or ( (isinstance(inchar,str)) and (inchar >= '0') and (inchar <= '9'))):
                        tmp = ' '
                        while (isinstance(inchar,str) and (inchar >= '0') and (inchar <= '9')):
                            tmp = tmp + inchar
                            inchar=fid.read(1)
                            #print inchar
                        if inchar == ',':
                            inchar =fid.read(1)
                            #print inchar
                        ndim = ndim +1
                        #print ndim
                        dim[ndim - 1] = int(tmp)
                        
                        
            elif (((inchar.isdigit()==False ) and (inchar != ',') and (inchar != '&') and (inchar != ' ') and (inchar != ':'))or (inchar=='') or (inchar=='_')):

                #rosacount += 1            
                #print rosacount               
                
                numitems=numitems+1
                nchar=0
#            # add in file.header...
                file[numitems-1,1]=int(offset)
                file[numitems-1,2]=int(ndim) 
#            file[numitems-1,6]=dim[:]    
                file[numitems-1,3]=int(dim[0])    
                file[numitems-1,4]=int(dim[1])    
                file[numitems-1,5]=int(dim[2])    
                



                newoff=1
                for i in range (0,ndim):
                    newoff=newoff*dim[i]
                offset=offset+newoff

                while (((inchar.isdigit()==False) and (inchar != ',') and (inchar != '&') and (inchar != ' ')and (inchar != ':')) or ((isinstance(inchar,str)) and (inchar>='0') and (inchar <= '9')) or (inchar==' ') or (inchar=='_')):
                    nchar=nchar+1
                    if (nchar<=200):    
                        tmpstr[nchar-1]=inchar
                    else:
                        return
                    inchar=fid.read(1)
                file[numitems-1,0]=(''.join(tmpstr[0:nchar])).lower()
                if (inchar==','):
                    inchar=fid.read(1)
            else:
                msg='Error reading header '+str(inchar)
                print msg

        inchar=fid.read(1)
        
        #print inchar
        
        while (inchar==' '):
            inchar =fid.read(1)
            #print inchar
#    
        fpos=(fid.tell() -1)
        
        #print fpos
#
        for i in range (0,numitems):
            file[i,1]=file[i,1]*4 + fpos
            
        #print numitems

#    return numitems

        #print file         
   
        return file
 
    
###################################
# END OF OPENBLASIUSFILE
###################################


def findBlasiusdata(file, dname):
#    Subroutine for readBlasiusdata.py whcih finds the variable loction
# within the file header.  


    import numpy as np
    
    lname = dname.lower()   #converts dname string to a string with all lower-case letters

    i=0
    ihdr=0
    #print numitems

    while (i<=numitems-1):
        if (dname == file[i,0]):
            ihdr=i
        i=i+1    
    
    ihdr=int(ihdr)
    return (ihdr)
    
###################################
# END OF FINDBLASIUSDATA
###################################

def readBlasiusdata(file, dname):

# Read in a varaible (dname) from FILE previously opened with
# openBlasiusFile.py. 
    
    import numpy as np

    ihdr = findBlasiusdata(file, dname)
    
#    if (ihdr ==0):
#        print 'Variable', dname, 'not found in header'
#        data=0
#        return

    fid.seek(file[ihdr,1])

    if (file[ihdr, 2]>=1):
        iip=file[ihdr,3]
#        print iip
    else:
        iip=1
    
    if (file[ihdr, 2]>=2):
        jjp=file[ihdr,4]
#        print jjp
    else:
        jjp=1
    
    if (file[ihdr, 2]>=3):
        kkp=file[ihdr,5]
#        print kkp
    else:
        kkp=1
    


    data=np.fromfile(fid,dtype='f', count=(iip*jjp*kkp), sep="")

    data=data.reshape((iip,jjp,kkp), order='F')        
            
    return data
    
    #print data


###################################
# END OF READBLASIUSDATA                        
###################################

def readgrids2d(file):

#  create 2d arrays for the coordinate grids.


    import numpy as np    

#    iip=readBlasiusdata(file, 'iip')
#    print iip
#    iip=int(iip)
#    jjp=readBlasiusdata(file, 'jjp')
#    jjp=int(jjp)
#    kkp=readBlasiusdata(file, 'kkp')
#    kkp=int(kkp)
    
                
    x=readBlasiusdata(file, 'x')
    xn=readBlasiusdata(file, 'xn')
    y=readBlasiusdata(file, 'y')
    yn=readBlasiusdata(file, 'yn')
    z=readBlasiusdata(file, 'z')
    zn=readBlasiusdata(file, 'zn')
    
    
    zsx=readBlasiusdata(file, 'zsx')
    zsy=readBlasiusdata(file, 'zsy')
    zsn=readBlasiusdata(file, 'zsn')
    
    zsnav=np.mean(zsn,axis=2)

    [xn2d,zn2d]=np.meshgrid(xn,zn)
    [x2d,z2d]=np.meshgrid(x,z)

    xn2d=xn2d.T
    zn2d=zn2d.T
    x2d=x2d.T
    z2d=z2d.T


# 

#    H=z[kkp-2,0]

#    for k in range (0,kkp-1):
#        zn2d[:,k]=zn2d[:,k]*(1-zsnav/H)+zsnav
#        z2d[:,k]=z2d[:,k]*(1-zsnav/H)+zsnav

    return xn2d, zn2d, x2d, z2d

########################################
#END OF READGRIDS2D
########################################


def levels(output_file):

        output_file_flip = output_file.T

        H = 1500

        Z1 = 0.25  

        R1 = 1.05 
    
        a0 = -Z1 
    
        a1 = 0.00
    
        a2 = Z1

        constant_spacing = 0.00

        Nlevels = len(output_file_flip)*2
        
        print Nlevels 

        levels = np.zeros(Nlevels)
        levels[0]=a0
        levels[1]=a1
        levels[2]=a2


        for i in range (3, Nlevels):
                # first get the required levels for the calculation
                a_minus_1 = levels[i-1]
                a_minus_2 = levels[i-2]
                a_n = 0    # initialising a_n, which gets calculated in the loop
    
                # test that current spacing multiplied by number of levels remaining is less than H
                remaining_levels = Nlevels - i
                remaining_height = H - a_minus_1
                constant_spacing = a_minus_1 - a_minus_2

                if (remaining_levels*constant_spacing <= remaining_height):
                #if true, then calculate next increment of height based on equation
                        a_n = a_minus_1 + (R1*(a_minus_1 - a_minus_2))
    
                else:    
                # otherwise use constant spacing
                        a_n = a_minus_1 + constant_spacing 

                        print "CONSTANT SPACING"
        
                        if a_n > H:
                                a_n = H
                # Now update the levels list to store current height, a_n
                levels[i]=a_n
    
        # Now need to extract Z levels (even) and ZN levels (odd)
        Z_levels=np.zeros(Nlevels/2)

        ZN_levels=np.zeros(Nlevels/2)

        for i in range(0,Nlevels/2):

                Z_levels[i] = levels[2*i+1]

                ZN_levels[i] = levels[2*i]

                print "Level " + str(i+1) + ", Height = " + str(Z_levels[i])

        return Nlevels, Z_levels

#-------------------------------------------------------------------------------#

def plot_data(u_data, w_data, p_data, Nlevels, Z_levels):
    
    
    flip_u = u_data.T   #flip the x and y axes
    
    flip_w = w_data.T

    flip_p = p_data.T

    #nflip = len(flip_u)

    #print nflip
    
    shorten_u =flip_u[1:80]   #remove first row as this is just zeros
    
    shorten_w =flip_w[1:80]
    
    shorten_p=flip_p[1:80]

    #nlevels = len(flip_u) #length of dataset
    
    #print nlevels
    
    #levels = list(range(0, nlevels)) #create list of numbers from 0 to the length of the data set (i.e. 0 to 49)
    
    #print levels
    
    mean_u = u_data.mean(axis=0)  

    lenMean = len(mean_u)

    print lenMean
    
    mean_w = w_data.mean(axis=0) 
    
    mean_p = p_data.mean(axis=0)

    #mean_u = shorten_u.mean(axis=1)   #calculate the mean wind speed for each level
    
    #mean_w = shorten_w.mean(axis=1)
   
    #mean_p = shorten_p.mean(axis=1)

    #print mean_w
    
    distance = 70    #choose the timeslice of speeds that you would like to compare against the mean
    
    slice = distance-1  
    
    #slice_u = shorten_u[:,slice]     #select the column corresponding to that timeslice
    
    #slice_w = shorten_w[:,slice]  
    
    #slice_p = shorten_p[:,slice]

    slice_u = flip_u[:,slice]

    slice_w = flip_w[:,slice]

    slice_p = flip_p[:,slice]


    
    
    #Create figure
    figure = plt.figure(1, facecolor = "white", figsize = (10,10))

    # Set up some basic parameters for the plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 8
    rcParams['legend.numpoints'] = 1
    axis_size = rcParams['font.size']+2
    
    #Add first axis
    
    #ax1 = figure.add_subplot(3,2,1)

    ax1 = plt.subplot2grid((24,2), (0,0), rowspan = 6)

    cax1 = ax1.imshow(shorten_u, aspect ='auto', cmap = "BuGn", origin = "lower")

    ax1.set_ylabel('Level', fontsize = axis_size)

    ax1.set_xlabel('Distance / m', fontsize = axis_size)

    ax1.annotate('a', xy=(0.05,0.93),  xycoords='axes fraction', backgroundcolor='none', horizontalalignment='left',  verticalalignment='top',  fontsize=axis_size)

    #cbar = figure.colorbar(cax1, ax=ax1, shrink=0.4, use_gridspec=True)

    #cbar = figure.colorbar(cax1, ax=ax1, orientation = "horizontal")
    
    ax7 = plt.subplot2grid((24,2), (6,0), rowspan = 2)

    cbar = plt.colorbar(cax1, cax = ax7, orientation = "horizontal")

    cbar.set_label('Streamwise velocity / m s$^{-1}$', multialignment='center', fontsize = axis_size)
    
    
    #Add third axis
    
    #ax3 = figure.add_subplot(3,2,3)

    ax3 = plt.subplot2grid((24,2), (8,0), rowspan = 6)

    cax3 = ax3.imshow(shorten_w, aspect ='auto', cmap = "PRGn", origin = "lower")

    ax3.set_ylabel('Level', fontsize = axis_size)

    ax3.set_xlabel('Distance / m', fontsize = axis_size)

    ax3.annotate('c', xy=(0.05,0.93),  xycoords='axes fraction', backgroundcolor='none', horizontalalignment='left',  verticalalignment='top',  fontsize=axis_size)

    #cbar = figure.colorbar(cax3, ax=ax3, shrink=0.4, use_gridspec=True)

    ax8 = plt.subplot2grid((24,2), (14,0), rowspan=2)

    cbar = plt.colorbar(cax3, cax = ax8, orientation = "horizontal")

    cbar.set_label('Vertical velocity / m s$^{-1}$', multialignment='center', fontsize = axis_size)
    
    
    #Add second axis
    
    #ax2 = figure.add_subplot(3,2,2)
    
    ax2 = plt.subplot2grid((24,2), (0,1), rowspan = 8)

    ax2.set_ylabel('Log height / m', fontsize = axis_size)

    ax2.set_xlabel('Streamwise velocity / m s$^{-1}$', fontsize = axis_size)

    ax2.plot(mean_u, Z_levels, color = "blue", ls = "--", label = "Mean")

    ax2.plot(slice_u, Z_levels, color = "green", label = "At x=%(distance)s" % {"distance" : distance})

    ax2.set_yscale('log')

    ax2.set_xlim(np.min(mean_u), np.max(mean_u)*1.15)

    ax2.set_ylim(0, np.max(Z_levels))

    plt.legend(loc = 'upper center')

    ax2.annotate('b', xy=(0.05,0.93),  xycoords='axes fraction', backgroundcolor='none', horizontalalignment='left',  verticalalignment='top',  fontsize=axis_size)
    

    #Add fourth axis
 
    #ax4 = figure.add_subplot(3,2,4)

    ax4 = plt.subplot2grid((24,2), (8,1), rowspan = 8)

    ax4.set_ylabel('Log height / m', fontsize = axis_size)

    ax4. set_xlabel('Vertical velocity / m s$^{-1}$', fontsize = axis_size)

    ax4.plot(mean_w, Z_levels, color = "blue", ls = "--", label = "Mean")

    ax4.plot(slice_w, Z_levels, color = "green", label = "At x=%(distance)s" % {"distance" : distance})

    ax4.set_yscale('log')

    plt.legend(loc = 'upper center')

    ax4.annotate('d', xy=(0.05,0.93),  xycoords='axes fraction', backgroundcolor='none', horizontalalignment='left',  verticalalignment='top',  fontsize=axis_size)

    majorLocator = MultipleLocator(0.01)

    ax4.xaxis.set_major_locator(majorLocator)


    #Add fifth axis

    #ax5 = figure.add_subplot(3,2,5)

    ax5 = plt.subplot2grid((24,2), (16,0), rowspan = 6)

    cax5 = ax5.imshow(shorten_p, aspect = 'auto',  cmap = "PRGn", origin = "lower")

    ax5.set_ylabel('Level', fontsize = axis_size)

    ax5.set_xlabel('Distance / m', fontsize = axis_size)

    ax5.annotate('e', xy=(0.05,0.93), xycoords='axes fraction', backgroundcolor='none', horizontalalignment='left', verticalalignment='top', fontsize=axis_size)

    ax9 = plt.subplot2grid((24,2), (22,0), rowspan = 2)

    cbar = plt.colorbar(cax5, cax=ax9, orientation="horizontal")

    cbar.set_label('Pressure / Pa', fontsize=axis_size)

    #cbaxes5 = figure.add_axes([0.8,0.1,0.03,0.8])

    #cb = plt.colorbar(ax5, cax=cbaxes5)

    #cbar = figure.colorbar(cax5, ax=ax5, shrink=0.4, use_gridspec=True, orientation = "horizontal")

    #cbar.set_label('Pressure / Pa', fontsize=axis_size)

    #ax5ins = inset_axes(ax5, width = "50%", height = "5%", loc = 4)

    #scalarMap._A = []

    #cbar = plt.colorbar(scalarMap._A, cax=ax5, orientation = "horizontal")

    #cbar.set_label("Pressure / Pa", fontsize=axis_size)

    #cbar.ax.xaxis.set_ticks_position("top")

    #cbar.ax.xaxis.set_label_position("top")

    #cbar.ax.tick_params(labelsize = axis_size, labelcolor = "white", color = "white")

    #cbar.outline.set_color("white")



    #Add sixth axis

    #ax6 = figure.add_subplot(3,2,6)

    ax6 = plt.subplot2grid((24,2), (16,1), rowspan = 8)

    ax6.set_ylabel('Log height / m', fontsize = axis_size)

    ax6.set_xlabel('Pressure /Pa', fontsize = axis_size)

    ax6.plot(mean_p, Z_levels, color = "blue", ls = "--", label = "Mean")

    ax6.plot(slice_p, Z_levels, color = "green", label="At x=%(distance)s"%{"distance" : distance})

    ax6.set_yscale('log')

    plt.legend(loc="upper center")

    ax6.annotate('f', xy=(0.05,0.93), xycoords='axes fraction', backgroundcolor='none', horizontalalignment='left', verticalalignment='top', fontsize=axis_size)


    
    #Make the plot look nicer
    
    plt.tight_layout(pad = 6, w_pad = 2, h_pad = 1)

    #Add a figure title

    figure.suptitle("2-d with canopy \n wake depth = 5", fontsize = 14)

    
    
    
    

if __name__ == "__main__":

    fname='/home/s1359318/rls_datastore/PhD/BLASIUS/blasius_v4.01/2-d-canopy/Output/RUN101_1025.DAT'
    
    u = 'u'
    
    v = 'v'
    
    w = 'w'
 
    p = 'p'
   
    out_file = openBlasiusfile(fname)
    
#    print out_file[:,0]  #print all the rows for the first column
    
    u_data = readBlasiusdata(out_file, u)[1:-1,1,:]
    
    v_data = readBlasiusdata(out_file, v)[1:-1,1,:]
    
    w_data = readBlasiusdata(out_file, w)[1:-1,1,:]
    
    p_data = readBlasiusdata(out_file, p)[1:-1,1,:]

    Nlevels, Z_levels = levels(u_data)

    plot_data(u_data, w_data, p_data, Nlevels, Z_levels)
    
   # plt.savefig('/home/s1359318/rls_datastore/PhD/BLASIUS/python/test/canopy/2dCNPY_wake5.png')
    
    plt.show()
    
    
    
