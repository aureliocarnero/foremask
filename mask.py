import os
import sys


from numpy import *
from  pylab import *

import matplotlib.pyplot as plt
import matplotlib.colors as col
matplotlib.rc('lines', linewidth=1.0)

import debug
import graphmask
from plots import scatter_plot


from scipy import array

matplotlib.rcParams['lines.linewidth'] = 3
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.labelsize'] = 'large'
matplotlib.rcParams['ytick.labelsize'] = 'large'
matplotlib.rcParams['figure.figsize'] = (10,8)
matplotlib.rcParams['figure.dpi'] = 40
matplotlib.rcParams['figure.autolayout']=True

def discrete_cmap(N=3):

#    """create a colormap with N (N<15) discrete colors and register it"""

    # define individual colors as hex values
    cpool = ['#000000', '#b20000', '#7FFF00']
    cmap3 = col.ListedColormap(cpool[0:N], 'indexed')
    return cmap3

def sanctify_gals(aninput,anoutput,io,logger):
#remove header line
    command='sed -i -e "1d" '+str(aninput)
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))

#Remove objects outside the mask (third column = blank)
    name="awk 'FS="
    name2='" " {if ($3!="") {print'
    name3=" $1,$2,$3;}else{print $1,$2,999;}}'"
    nametot=str(name)+str(name2)+str(name3)+' '+str(aninput)+' > '+str(anoutput)
    os.system(nametot)
    return anoutput


def prephealpixtile(res,trimmasks,bolemask,ramin,ramax,decmin,decmax,io,logger):
    import healpix as hpp
    name = 'pixels_coordinates_%s'%str(res)
    hpp.healRDfullsky(res,name)
    ra,dec = loadtxt(name,unpack=True)
    heal_in = arange(len(ra))
    ma = (ra>=ramin)*(ra<=ramax)*(dec>=decmin)*(dec<=decmax)
    ra=ra[ma]
    dec=dec[ma]
    heal_in = heal_in[ma]
    savetxt(name,array([ra,dec]).transpose(),fmt='%e',delimiter='\t')
    healpreptile = []
    for i,mask in enumerate(trimmasks):
        if bolemask[i]:
            outname = str(mask)+'_healpix'
            command='polyid -W %s %s %s'%(mask,name,'temp.in')
            mvcount=io.run_command(command,logger)
            logger.info("Executed %s. Return code = %s" % (command,mvcount))
            sanctify_gals('temp.in',outname,io,logger)
            ra,dec,w = loadtxt(outname,unpack=True)
            ma = (w!=999)
            heal_i = heal_in[ma]
            w = w[ma]
            savetxt(outname,array([heal_i,w]).transpose(),fmt='%i %10.6f',delimiter='\t')
            healpreptile.append(outname)
        else:
            healpreptile.append('nomask')
    return healpreptile


def paintmask_fromlist_redgreen(lisst,maskout,title,io,logger):
    cmap3 = discrete_cmap()
    bounds =[0,1,2,3]
    norm = col.BoundaryNorm(bounds, cmap3.N)
    #mvcount=io.run_command(command,logger)
    
    #logger.info("Executed %s. Return code = %s" % (command,mvcount))
    
    cp,m = graphmask.plot_mangle_map(lisst,autoscale=True,enlarge_border=0.01,edgecolor='black',draw_colorbar=False,plottitle=str(title),bgcolor='white',cmap=cmap3,vmin=0.,vmax=3.,norm=norm)
    cbar = colorbar(cp, ticks=[0.5,1.5,2.5])
    cbar.ax.set_yticklabels(['error', 'failed', 'passed'])

    
    
    plt.savefig(maskout)
    
    plt.close()

def create_mangle_map(io,logger,ramin,ramax,decmin,decmax,idds,ww):
    f = open('poly.pol','w')

    f.write('rectangle\n')
    f.write('unit d\n')
    for mi,ma,dmi,dma,w in zip(ramin,ramax,decmin,decmax,ww):
        
        
        f.write('%f %f %f %f %f\n'%(mi,ma,dmi,dma,w))
        
    f.close()
    
    command='poly2poly -ir1d -ol poly.pol poly.list'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))
    aa,bb,cc,dd = loadtxt('poly.list.weight',unpack=True)
    f = open('poly.list.weight','w')
    for a,w,c,d in zip(aa,ww,cc,dd):
        f.write('%d %f %f %f\n'%(a,w,c,d))
    f.close()
    return 'poly.list','poly.list.weight'

def preptimetile(filter_list,bole,tilename,mask_path,io,logger):
   

    maskplots = []
    masktime = []
    www = []
    for i,filter in enumerate(filter_list):
        if bole[i]:
            print bole[i]
            maskmask = str(mask_path)+str(tilename)+'_holymolys_time_'+str(filter)+'.pol'
            maskout = 'trimmask_time_'+str(filter)+'.pol'
            command='cp '+str(maskmask)+' '+str(maskout)
            mvcount=io.run_command(command,logger)
            logger.info("Executed %s. Return code = %s" % (command,mvcount))
            masktime.append(maskout)
	    listtemp = str(mask_path)+str(tilename)+'_holymolys_time_'+str(filter)+'.list'

            listout = 'trimmask_time_'+str(filter)+'.list'

	    if os.path.isfile(listtemp): 
		command = 'cp '+listtemp+' '+listout
                mvcount=io.run_command(command,logger)
                logger.info("Executed %s. Return code = %s" % (command,mvcount))
		command = 'cp '+listtemp+'.weight '+listout+'.weight'

	    else:
                command='poly2poly -ol30 '+str(maskout)+' '+str(listout)
            mvcount=io.run_command(command,logger)
            logger.info("Executed %s. Return code = %s" % (command,mvcount))

            cp,m = graphmask.plot_mangle_map(listout,autoscale=True,enlarge_border=0.01,bgcolor='white',cmap=cm.jet)
            plotname= 'trimmask_time_'+str(filter)+'.png'
            plt.savefig(plotname)

            plt.close()
            maskplots.append(plotname)
            os.system('grep weight '+str(maskout)+' > tt')

            name="awk 'FS="
            name2='" " {print'
            name3=" $6;}'"
            nametot=str(name)+str(name2)+str(name3)+' tt > time_'+str(filter)+'.w'
            os.system(nametot)
            www.append('time_'+str(filter)+'.w')
        else:
            maskplots.append("nomask")
            masktime.append("nomask")
            www.append("nomask")


    return maskplots,masktime,www



def preparationtile(filter_list,tilename,mask_path,io,logger):
    
   
    maskmaglim = []
    masklist = []
    maskplots = []
    masklistweight = []
    www = []
    boolmask = []
    for filter in filter_list:
        maskmask = str(mask_path)+str(tilename)+'_holymolys_maglims_'+str(filter)+'.pol'
        maskout = 'trimmask_maglim_'+str(filter)+'.pol'
#The actual triming
        try:
            command='cp '+str(maskmask)+' '+str(maskout)
            mvcount=io.run_command(command,logger)
            logger.info("Executed %s. Return code = %s" % (command,mvcount))	
            maskmaglim.append(maskout)
	    listtemp = str(mask_path)+str(tilename)+'_holymolys_maglims_'+str(filter)+'.list'

            listout = 'trimmask_maglim_'+str(filter)+'.list'
	    if os.path.isfile(listtemp):
		command='cp '+listtemp+' '+listout
		mvcount=io.run_command(command,logger)
		logger.info("Executed %s. Return code = %s" % (command,mvcount))
 		command='cp '+listtemp+'.weight '+listout+'.weight'
		mvcount=io.run_command(command,logger)
		logger.info("Executed %s. Return code = %s" % (command,mvcount))

	    else:
#create list and paint the mask
                command='poly2poly -ol30 '+str(maskout)+' '+str(listout)
                mvcount=io.run_command(command,logger)
                logger.info("Executed %s. Return code = %s" % (command,mvcount))

            masklist.append(listout)
            listweight = str(listout)+'.weight'
            masklistweight.append(listweight)
            
            
            #cp,m = graphmask.plot_mangle_map(listout,autoscale=True,enlarge_border=0.01,bgcolor='white',cmap=cm.jet)

            plotname= 'trimmask_maglim_'+str(filter)+'.png'

            command = 'cp '+str(mask_path)+str(tilename)+'_holymolys_maglims_'+str(filter)+'.png '+plotname
            mvcount=io.run_command(command,logger)
            logger.info("Executed %s. Return code = %s" % (command,mvcount))

    
            #plt.savefig(plotname)
    
            #plt.close()
            maskplots.append(plotname)
            os.system('grep weight '+str(maskout)+' > tt')

            name="awk 'FS="
            name2='" " {print'
            name3=" $6;}'"
            nametot=str(name)+str(name2)+str(name3)+' tt > maglim_'+str(filter)+'.w'
            os.system(nametot)
            www.append('maglim_'+str(filter)+'.w')
            boolmask.append(True)
        except:
            boolmask.append(False)
            maskmaglim.append("nomask")
            masklist.append("nomask")
            masklistweight.append("nomask")
            maskplots.append("nomask")
            www.append("nomask")




    return maskmaglim,masklist,masklistweight,maskplots,www,boolmask


def getstatsflags_preparation(trimmasks,galaxycat,io,logger): 
    
    
    wie = []
    num_in = []
    num_out = []
    num_masked = []
    mak_in = []
    
    for maskis in trimmasks:
        if maskis != "nomask":    
            outcatflag = 'temp_'+str(maskis)
            print maskis,outcatflag
#match galaxies in galaxy catalog with polygon weight from the given mask
            command='polyid -W '+str(maskis)+' '+str(galaxycat)+' '+str(outcatflag)
            mvcount=io.run_command(command,logger)
            logger.info("Executed %s. Return code = %s" % (command,mvcount))	


#remove header line and objects outside the mask
            sanctify_gals(outcatflag,'gal.in',io,logger)


            r,d,nu=loadtxt('gal.in',unpack=True)
            wie.append(nu)
            os.remove(str(outcatflag))
            os.remove('gal.in')
        else:
            ra,dec = loadtxt(galaxycat,unpack=True)
            xxx = len(ra)
            wie.append([999]*xxx)


    totalmask = [True] * len(wie[0])
    
    for w,maskis in zip(wie,trimmasks):
        if maskis != "nomask": 
            makmask=(w==0)
            makin=(w!=0)*(w!=999)
            mak_in.append(makin)
            makout=(w==999)
            numin=len(r[makin])
        
            numout=len(r[makout])
            numask=len(r[makmask])

            num_in.append(numin)
            num_out.append(numout)
            num_masked.append(numask)
            totalmask = totalmask*makin

        else:
            mak_in.append([False]*len(r))
            num_in.append(0)
            num_out.append(len(r))
            num_masked.append(0)
#    totalmask = [False]*len(r)
            
    print len(r[totalmask])

#Return number of galaxies inside, outside and masked out
    os.remove(galaxycat)
    try:
        return len(r),len(r[totalmask]),totalmask,num_in,num_out,num_masked
    except:
        return len(r),len(r),totalmask,num_in,num_out,num_masked



def run_properties_magarea(mask,band,io,logger):


    
    


    os.system('grep weight '+str(mask)+' > tt')

    name="awk 'FS="
    name2='" " {print'
    name3=" $6;}'"
    nametot=str(name)+str(name2)+str(name3)+' tt > maglims_'+str(band)

    

    os.system(nametot)            

    name3=" $10;}'"

    nametot=str(name)+str(name2)+str(name3)+' tt > area_'+str(band)

    os.system(nametot)

   




    name2,totarea,effearea,blankarea,mag95=areavsmag('maglims_'+str(band),'area_'+str(band),band)
    return name2,totarea,effearea,blankarea,mag95

def densityvsmag(mask_name,mask_out,X,io,logger):

    
    command = 'grep weight '+str(mask_name)+' > tempa'
    os.system(command)

    command1 = "awk 'FS="
    command1 += '" " {print $2;}'
    command1 +="' tempa > polygons"
    os.system(command1)

    command1 = "awk 'FS="
    command1 += '" " {print $6;}'
    command1 +="' tempa > maglims"
    os.system(command1)



    polys = loadtxt("polygons",unpack=True)
    magiis = loadtxt("maglims",unpack=True)
    n=len(polys)

    
    ra = X.data.field('ra')
    dec = X.data.field('dec')
    savetxt("cat.in",array([ra,dec]).transpose(),fmt='%e',delimiter='\t')

 
    command='polyid '+str(mask_name)+' cat.in tempid'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))

    sanctify_gals('tempid','galid.in',io,logger)


    getareamask(mask_name,'area.w',io)

    areas = loadtxt("area.w",unpack=True)

    idgal = loadtxt("galid.in",unpack=True)
    numobjs = []
    dfile = open("number_.w", "w")
    for poll in polys:

        numobjs.append(idgal[idgal == poll].size)
        dfile.write(str(idgal[idgal == poll].size)+'\n')
    dfile.close()
    dfile = open("density_.w", "w")
 
    for i,val in enumerate(areas):

        #if ((areas[i] != 0)& (areas[i] < 1.)):
        if ((numobjs[i]>0)&(magiis[i]>0)):

            #densit = log(float(numobjs[i])/(float(areas[i])*1.18181034378e+07))
            densit = float(numobjs[i])/(float(areas[i])*1.18181034378e+07)
 
        else:
            densit = 0.0
        #else:
        #    densit = 0.0
        #if (densit<1000.):
        dfile.write(str(densit)+'\n')
        #else:
            #dfile.write(str(0.0)+'\n')

    dfile.close()
    command = 'weight -zdensity_.w '+str(mask_name)+' '+str(mask_out)
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))

    #return "density_"+str(band)+".w"

#Remove objects blanked by mask, weight=0. Substitute weight of all objects to 1 
    

    


#no preciso hacer lo de las galaxias con area, tengo que coger el area como un peso de cada polygon en orden.
#Despues, histograma de id y convertirlo en un peso para cada polygon, dividir uno por otro y ya tengo el peso para la mascara. Coger la trim mask ya! quizas haya que modificar alguna cosa. va a ser foda!
#aceptar molly y nacho como usuarios para ver, junto a DESDM
#also in product log, information about the mask origin and date





def magvsexp(explist,maglist,band,io,logger):
    exposures=loadtxt(explist,unpack=True)
    maglims=loadtxt(maglist,unpack=True)
    mas=(maglims!=0)
    pname='exposurevsmag_'+str(band)+'.png'
    scatter_plot(pname,exposures[mas],maglims[mas],'Exposure time','Magnitude limit',axislim=False)
    return pname

def maskmaglim4ran(maskname,maglim,maskout,io,logger): 

    command1 = 'grep weight '+str(maskname)+' | awk '
    command2 = "'FS="
    command3 = '" " {print $6}'
    command4 = "' > maglist"
    command = str(command1)+str(command2)+str(command3)+str(command4)
    os.system(command)
    command1 = "awk '{if ($1<"+str(maglim)+"){print 0;}else{print 1}}' maglist > weightmaglim"
    os.system(command1)
    command = 'weight -zweightmaglim '+str(maskname)+' '+str(maskout)
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))


def getstats(maskname,X,io,logger): 

    ra = X.data.field('ra')
    dec = X.data.field('dec')
    galaxycat='makestats'
    savetxt(galaxycat,array([ra,dec]).transpose(),fmt='%e',delimiter='\t')




#match galaxies in galaxy catalog with polygon weight from the given mask
    command='polyid -W '+str(maskname)+' '+str(galaxycat)+' temp'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))	

    sanctify_gals('temp','gal.in',io,logger)


    r,d,nu=loadtxt('gal.in',unpack=True)
    makmask=(nu==0)
    makin=(nu!=0)*(nu!=999)
    makout=(nu==999)
    numin=len(r[makin])
    numout=len(r[makout])
    numask=len(r[makmask])

#Remove temporary files
    command='rm temp'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))

#    command='rm '+str(galaxycat)
#    mvcount=io.run_command(command,logger)
#    logger.info("Executed %s. Return code = %s" % (command,mvcount))

    command='rm gal.in'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))


    command='rm '+str(galaxycat)
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))

#Return number of galaxies inside, outside and masked out

    return numin,numout,numask

#THIS MODULE PAINT MAGLIM VERSUS AREA FOR A GIVEN BAND AND A GIVEN MASK
def areavsmag(maglim,magarea,band):
    try:
        area=loadtxt(magarea,unpack=True)
        mag=loadtxt(maglim,unpack=True)
        mak=(mag!=0)



        totalarea=sum(area)*3283
        effearea=sum(area[mak])*3283
        blankarea=(totalarea-effearea)


        minmag=min(mag[mak])
        maxmag=max(mag[mak])
        pos=linspace(minmag,maxmag,100)
        aaa = [effearea]*100
                                                   
                                                       
                                                           #events, edges, patches = plt.hist(mag[mak],bins=histt,histtype='stepfilled')
                                                               #print events,edges,patches
        makmag=mag[mak]
        makarea=area[mak]
       
        d=digitize(makmag,pos)
        for i,vals in enumerate(makmag):
            indd=d[i]

            aaa[indd:100]=aaa[indd:100]-(makarea[i]*3283)
                                                                                                                   

        width = (maxmag-minmag)/100.
        
        inde = min(range(len(aaa)), key=lambda i: abs(aaa[i]-(0.90*effearea)))

        figure(4)#,figsize=(6.4,5))
        ax1 = subplot(111)

        plt.bar(pos,aaa,width,color='b')
        plt.axhline(y=(0.9*effearea), linewidth=1,color='r')
        plt.axvline(x=pos[inde], linewidth=1,color='r')
        plt.ylabel('Effective Area [square degrees]')
        plt.xlabel('Magnitude Limit')
        plt.ylim((0,totalarea))
        ax1.xaxis.tick_bottom()
        ax1.yaxis.tick_left()
        ax1.text((0.99*maxmag),(0.91*effearea),'90 %',fontsize=12, color='red')

        
        

        
        ax2 = twinx()
        plt.ylabel('Effective Area [%]')
        ax2.yaxis.tick_right()
        plt.ylim((0,(100*totalarea/effearea)))
        

        namefig='area_vs_mag_'+str(band)+'.png'
        plt.savefig(namefig)
        close(4)
        return namefig,totalarea,effearea,blankarea,pos[inde]
    except:
        namefig='area_vs_mag_'+str(band)+'.png'
        plt.savefig(namefig)
        close()
        totalarea = 0
        effearea = 0
        blankarea = 0
        
        return namefig,totalarea,effearea,blankarea,0

  

#THIS MODULE TRIM A WHOLE MASK WITHIN A CERTAIN RA-DEC RANGE
def trim_mask(maskmask,ramin,ramax,decmin,decmax,namm,io,logger):


#Create a polygon file defining the region where to trim the mask and pixelize to the DES mask resolution
    dfile = open("region.rect", "w")
    newLine ='rectangle \n unit d \n'+str(ramin)+' '+str(ramax)+' '+str(decmin)+' '+str(decmax)+'\n'  
    dfile.write(newLine)
    dfile.close()

    command='pixelize -Ps0,9 region.rect region.pol'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))	

    command='snap region.pol region_snapped.pol'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))	





#The actual triming
    command='trim_mask.sh '+str(maskmask)+' region_snapped.pol '+str(namm)
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))	

#Delete temporary files	
    command='rm region.pol'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))


    command='rm region_snapped.pol'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))
def trim_mask_component(maskmask,ramin,ramax,decmin,decmax,namm):


#Create a polygon file defining the region where to trim the mask and pixelize to the DES mask resolution
    dfile = open("region.rect", "w")
    newLine ='rectangle \n unit d \n'+str(ramin)+' '+str(ramax)+' '+str(decmin)+' '+str(decmax)+'\n'  
    dfile.write(newLine)
    dfile.close()

    command='pixelize -Ps0,4 region.rect region.pol'
    os.system(command)

#    mvcount=io.run_command(command,logger)
#    logger.info("Executed %s. Return code = %s" % (command,mvcount))	

    command='snap region.pol region_snapped.pol'
    os.system(command)

#    mvcount=io.run_command(command,logger)
#    logger.info("Executed %s. Return code = %s" % (command,mvcount))	

#The actual triming
    command='trim_mask.sh '+str(maskmask)+' region_snapped.pol '+str(namm)
    os.system(command)


def trim_mask_list(maskmask,region,io,logger):

    command='pixelize -Ps0,9 '+str(region)+' region.pol'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))

    command='snap region.pol region_snapped.pol'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))

    command='trim_mask.sh '+str(maskmask)+' region_snapped.pol trimregion.pol'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))

    command='rm region.pol'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))

    command='rm region_snapped.pol'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))



#THIS MODULE SELECT GALAXIES THROUGH MASK. ONCE THE MASK DATABASE IS CREATED IN THE INSTALLATION PROCEDURE, THIS MODULE WILL BE OBSOLETE
def select_galaxies(maskname,galaxycat,galaxyout,io,logger): 





#match galaxies in galaxy catalog with polygon weight from the given mask
    command='polyid -W '+str(maskname)+' '+str(galaxycat)+' temp'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))	

    sanctify_gals('temp','gal.in',io,logger)


#Remove objects blanked by mask, weight=0. Substitute weight of all objects to 1 
    name="awk 'FS="
    name2='" " {if ($3!=0) {print'
    name3=" $1,$2,1;}}'"
    nametot=str(name)+str(name2)+str(name3)+' gal.in > '+str(galaxyout)
    print str(nametot)
    os.system(nametot)

#Read number of galaxies
    ratemp,dectemp,ttt=loadtxt(galaxyout,unpack=True)
    n_gals=len(ratemp)

#Remove temporary files
    command='rm temp'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))

#    command='rm '+str(galaxycat)
#    mvcount=io.run_command(command,logger)
#    logger.info("Executed %s. Return code = %s" % (command,mvcount))

    command='rm gal.in'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))

#Return number of galaxies
    return n_gals

#CREATE RANDOM CATALOG BASED IN A MASK
def create_random(n_randoms,maskname,maskout,io,logger):    

#Create random catalog
    command='ransack -c0 -r'+str(n_randoms)+' '+str(maskname)+' trimregion.random'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))
#remove header line
    command='sed -i -e "1d" trimregion.random'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))

#create random catalog, with weight=1
    ran_ra,ran_dec,ran_id=loadtxt('trimregion.random',unpack=True)
    savetxt(maskout,array([ran_ra,ran_dec,ones(n_randoms)]).transpose(),fmt='%e',delimiter='\t')

#remove temporary files
    command='rm trimregion.random'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))


def paintt(mask,io,logger):
    
    command='poly2poly -ol30 '+str(mask)+' trimregion.list'
   
    mvcount=io.run_command(command,logger)
    
    logger.info("Executed %s. Return code = %s" % (command,mvcount))

    cp,m = graphmask.plot_mangle_map('trimregion.list',autoscale=True,enlarge_border=0.01,bgcolor='white',cmap=cm.jet)

    name='mask.png'
    
    plt.savefig(name)
    
    plt.close()

    return name
def paintmask(mask,maskout,io,logger):
    
    command='poly2poly -ol30 '+str(mask)+' trimregion.list'
    os.system(command)   
    #mvcount=io.run_command(command,logger)
    
    #logger.info("Executed %s. Return code = %s" % (command,mvcount))

    cp,m = graphmask.plot_mangle_map('trimregion.list',autoscale=True,enlarge_border=0.01,bgcolor='white',cmap=cm.jet)

    
    
    plt.savefig(maskout)
    
    plt.close()

def paintmask_consolida(mask,rm,rM,dm,dM,maskout,io,logger):
    
    command='poly2poly -ol30 '+str(mask)+' trimregion.list'
    os.system(command)   
    #mvcount=io.run_command(command,logger)
    
    #logger.info("Executed %s. Return code = %s" % (command,mvcount))

    cp,m = graphmask.plot_mangle_map('trimregion.list',autoscale=True,enlarge_border=0.01,minaz=rm,maxaz=rM,minel=dm,maxel=dM,bgcolor='white',cmap=cm.jet)

    
    
    plt.savefig(maskout)
    
    plt.close()


def paintmasklist(list,plotout,io,logger):
    

    cp,m = graphmask.plot_mangle_map(list,autoscale=True,enlarge_border=0.01,bgcolor='white',cmap=cm.jet)

    
    
    plt.savefig(plotout)
    
    plt.close()
def paintmasklistwithpoints(list,ra,dec,plotout,io,logger):
    
    plt.clf()
    plt.figure(1)#, figsize=(8,6))
    plt.plot(ra,dec,'ko',alpha=0.75)
          


    cp,m = graphmask.plot_mangle_map(list,autoscale=True,enlarge_border=0.01,bgcolor='white',cmap=cm.jet)

    
    
    plt.savefig(plotout)
    
    plt.close(1)

def paintmaskdens(mask,plotout,io,logger):
    
    command='poly2poly -ol30 '+str(mask)+' trimregion.list'
    #os.system(command)   
    mvcount=io.run_command(command,logger)
    
    logger.info("Executed %s. Return code = %s" % (command,mvcount))
    
    cp,m = graphmask.plot_mangle_map('trimregion.list',vmax=8,autoscale=True,enlarge_border=0.01,bgcolor='white',cmap=cm.winter)

    
    
    plt.savefig(plotout)
    
    plt.close()
def getareamask(maskin,weightout,io):
       

    os.system('grep weight '+str(maskin)+' > tt')

    name="awk 'FS="
    name2='" " {print'
    name3=" $10;}'"
    nametot=str(name)+str(name2)+str(name3)+' tt > '+str(weightout)

    

    os.system(nametot)            
    

   
def getweightsmask(maskin,weightout,io):
    print 'empieza getweightsmask'    


    

    os.system('grep weight '+str(maskin)+' > tt')

    name="awk 'FS="
    name2='" " {print'
    name3=" $6;}'"
    nametot=str(name)+str(name2)+str(name3)+' tt > '+str(weightout)

    

    os.system(nametot)            
    print 'termina getweightsmask con',weightout

def createweightsconsolidator(band,release,io,logger):
    


    said='trimregion_'+str(band)+'_maglims.pol'

    os.system('grep weight '+str(said)+' > tt')

    name="awk 'FS="
    name2='" " {print'
    name3=" $6;}'"
    nametot=str(name)+str(name2)+str(name3)+' tt > maglims_'+str(band)

    

    os.system(nametot)            
    mask='/archive/staging/DES/coadds/masks/'+str(release)+'/'+str(release)+'_holymolys_time_'+str(band)+'.pol'
    
    said='trimregion_'+str(band)+'_time.pol'
    trim_mask_list(mask,'region.rect',io,logger)
	
    command = 'mv trimregion.pol '+str(said)
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))


    os.system('grep weight '+str(said)+' > tt')

    
    nametot=str(name)+str(name2)+str(name3)+' tt > time_'+str(band)

    os.system(nametot)

def createweights(band,release,io):
    logger = lib.log.get_logger()
#    if io == None:
#        io, debuglog = debug.beginLog('Mask')
 
    
    
    try:
        ramin=io.getConfigById("ramin")
        ramax=io.getConfigById("ramax")
        decmin=io.getConfigById("decmin")
        decmax=io.getConfigById("decmax")   
    except Exception,e:
        print e
        ramin=io.getConfigById("RAMIN")
        ramax=io.getConfigById("RAMAX")
        decmin=io.getConfigById("DECMIN")
        decmax=io.getConfigById("DECMAX")   
    mask='/archive/staging/DES/coadds/masks/'+str(release)+'/'+str(release)+'_holymolys_maglims_'+str(band)+'.pol'
    said='trimregion_'+str(band)+'_maglims.pol'
    trim_mask(mask,ramin,ramax,decmin,decmax,said,io,logger)
    os.system('grep weight '+str(said)+' > tt')
    #comfour = 'grep weight trimregion.pol > tt'
    #mvcount=io.run_command(comfour,logger)
    #io.logger.info("Executed %s. Return code = %s" % (comfour,mvcount))
    name="awk 'FS="
    name2='" " {print'
    name3=" $6;}'"
    nametot=str(name)+str(name2)+str(name3)+' tt > maglims_'+str(band)
    #mvcount=io.run_command(nametot,stdout='maglims')
    #io.logger.info("Executed %s. Return code = %s" % (comfour,mvcount))
    

    os.system(nametot)            
    mask='/archive/staging/DES/coadds/masks/'+str(release)+'/'+str(release)+'_holymolys_time_'+str(band)+'.pol'
    
    said='trimregion_'+str(band)+'_time.pol'
    trim_mask(mask,ramin,ramax,decmin,decmax,said,io,logger)
	



    os.system('grep weight '+str(said)+' > tt')

    
    nametot=str(name)+str(name2)+str(name3)+' tt > time_'+str(band)

    os.system(nametot)


def run(band,release,io):
#    logger = lib.log.get_logger()
#    if io == None:
#        io, debuglog = debug.beginLog('Mask')
    logger = lib.log.get_logger()


    
    
    try:
        ramin=io.getConfigById("ramin")
        ramax=io.getConfigById("ramax")
        decmin=io.getConfigById("decmin")
        decmax=io.getConfigById("decmax")   
    except Exception,e:
        print e
        ramin=io.getConfigById("RAMIN")
        ramax=io.getConfigById("RAMAX")
        decmin=io.getConfigById("DECMIN")
        decmax=io.getConfigById("DECMAX")   
    mask='/archive/staging/DES/coadds/masks/'+str(release)+'/'+str(release)+'_holymolys_maglims_'+str(band)+'.pol'
    said='trimregion_'+str(band)+'_maglims.pol'
    trim_mask(mask,ramin,ramax,decmin,decmax,said,io,logger)
    
    command='poly2poly -ol30 '+str(said)+' trimregion.list'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))



    cp,m = graphmask.plot_mangle_map('trimregion.list',autoscale=True,enlarge_border=0.01,bgcolor='white',cmap=cm.jet)
    nameI='mask_'+str(release)+'_holymolys_maglims_'+str(band)+'.png'
    plt.savefig(nameI)
    plt.close()

    os.system('grep weight '+str(said)+' > tt')
    name="awk 'FS="
    name2='" " {print'
    name3=" $6;}'"
    nametot=str(name)+str(name2)+str(name3)+' tt > maglims_'+str(band)
    #mvcount=io.run_command(nametot,stdout='maglims')
    #io.logger.info("Executed %s. Return code = %s" % (comfour,mvcount))
    

    os.system(nametot)            
                            
    
    name3=" $10;}'"
    #io.logger.info("Executed %s. Return code = %s" % (comfour,mvcount))
    nametot=str(name)+str(name2)+str(name3)+' tt > area_'+str(band)

    os.system(nametot)

    
    #name2='mask_'+str(release)+'_area_vs_maglims_'+str(band)+'.png'
    name2,totarea,effearea,blankarea=areavsmag('maglims_'+str(band),'area_'+str(band),band)
    return nameI,name2,totarea,effearea,blankarea
def runconsolidator(band,release,maskfortrim,io):
#    logger = lib.log.get_logger()
#    if io == None:
#        io, debuglog = debug.beginLog('Mask')
    logger = lib.log.get_logger()


    
    
    try:
        ramin=io.getConfigById("ramin")
        ramax=io.getConfigById("ramax")
        decmin=io.getConfigById("decmin")
        decmax=io.getConfigById("decmax")   
    except Exception,e:
        print e
        ramin=io.getConfigById("RAMIN")
        ramax=io.getConfigById("RAMAX")
        decmin=io.getConfigById("DECMIN")
        decmax=io.getConfigById("DECMAX")   
#    mask='/archive/staging/DES/coadds/masks/'+str(release)+'/'+str(release)+'_holymolys_maglims_'+str(band)+'.pol'
    said='trimregion_'+str(band)+'_maglims.pol'
#    trim_mask(mask,ramin,ramax,decmin,decmax,said,io,logger)
    
    command='poly2poly -ol30 '+str(said)+' trimregion.list'
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))


    #command='poly2poly -ol30 trimregion.pol trimregion.list'
    #mvcount=io.run_command(command,logger)
    #io.logger.info("Executed %s. Return code = %s" % (command,mvcount))

    cp,m = graphmask.plot_mangle_map('trimregion.list',autoscale=True,enlarge_border=0.01,bgcolor='white',cmap=cm.jet)
    nameI='mask_'+str(release)+'_holymolys_maglims_'+str(band)+'.png'
    plt.savefig(nameI)
    plt.close()

    os.system('grep weight '+str(said)+' > tt')
    #comfour = 'grep weight trimregion.pol > tt'
    #mvcount=io.run_command(comfour,logger)
    #io.logger.info("Executed %s. Return code = %s" % (comfour,mvcount))
    name="awk 'FS="
    name2='" " {print'
    name3=" $6;}'"
    nametot=str(name)+str(name2)+str(name3)+' tt > maglims_'+str(band)
    #mvcount=io.run_command(nametot,stdout='maglims')
    #io.logger.info("Executed %s. Return code = %s" % (comfour,mvcount))
    

    os.system(nametot)            
    #comfour = 'awk \'FS=\" \" {printi $6}\' tt > maglims'
    #mvcount=io.run_command(comfour,logger)
    #io.logger.info("Executed %s. Return code = %s" % (comfour,mvcount))
                            
    #comone = "grep weight trimregion.pol | awk \'FS=\" \" {print $6}\' > maglims"
    
    #mvcount=io.run_command(comone,logger)
    #io.logger.info("Executed %s. Return code = %s" % (comone,mvcount))
    mask='/archive/staging/DES/coadds/masks/'+str(release)+'/'+str(release)+'_molys_area_'+str(band)+'.pol'
    #mask='/archive/staging/DES/coadds/masks/'+str(release)+'mask/'+str(release)+'_molys_area_g.pol'
    said='trimregion_'+str(band)+'_area.pol'

    trim_mask_list(mask,maskfortrim,io,logger)

    comfour='mv trimregion.pol '+str(said)
    mvcount=io.run_command(comfour,logger)
    logger.info("Executed %s. Return code = %s" % (comfour,mvcount))


    os.system('grep weight '+str(said)+' > tt')

    #comfour='grep weight trimregion.pol > tt'
    #mvcount=io.run_command(comfour,logger)
    #io.logger.info("Executed %s. Return code = %s" % (comfour,mvcount))
    nametot=str(name)+str(name2)+str(name3)+' tt > area_'+str(band)

    os.system(nametot)

    #comfour="awk \'FS=\" \" {print $6}\' tt > area"
    #mvcount=io.run_command(comfour,logger)
    #io.logger.info("Executed %s. Return code = %s" % (comfour,mvcount))
    



    #name2='mask_'+str(release)+'_area_vs_maglims_'+str(band)+'.png'
    name2,totarea,effearea,blankarea=areavsmag('maglims_'+str(band),'area_'+str(band),band)
    return nameI,name2,totarea,effearea,blankarea
#    io.beginSection('test')
    #io.beginSubSection(str(band),indent_level='10')
    #io.addPlot(nameI, 'Magnitude limit', 'Magnitude limit mask in mag_aper2 at 10sigma in band ' + str(band)+' for '+str(release)+'.')
    #io.addPlot(name2, 'Magnitude limit vs Area', 'Effective area as a function of magnitude limit in mag_aper2 at 

def mergemask(listmask,maskout,idd,io,logger,name=None):
    weightss = []
    for i,maskis in enumerate(listmask):
        wei='weight_'+str(i)
        weightss.append(wei)
        getweightsmask(maskis,wei,io)
    if name:
        weei='totweights_'+str(name)+'_'+str(idd)
    else:
        weei='totweights_'+str(idd)
    command='cat '
    for ww in weightss:
        command +=str(ww)+' '
    command +='> '+str(weei)
    os.system(command)
    command = 'weight -z'+str(weei)+' '
    for i,maskis in enumerate(listmask):
        command +=str(maskis)+' '
    command +=str(maskout)
    mvcount=io.run_command(command,logger)
    logger.info("Executed %s. Return code = %s" % (command,mvcount))
    return weei




