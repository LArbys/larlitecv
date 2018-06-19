import sys,os,commands

dirs=[]
for d in os.listdir(os.environ['LARLITECV_BUILDDIR']):
    if not len([x for x in os.listdir('%s/%s' % (os.environ['LARLITECV_BUILDDIR'],d)) if x.endswith('.o')]): continue
    dirs.append(d)
libs=[x for x in commands.getoutput('larlitecv-config --libs').split() if not x.startswith('-llarlitecv')]
libs+= commands.getoutput('root-config --libs').split()
if 'LARLITE_BASEDIR' in os.environ:
    libs+= commands.getoutput('larlite-config --libs').split()
    libs+= ["-lBasicTool_FhiclLite","-lBasicTool_GeoAlgo","-lSelectionTool_OpT0FinderAna","-lSelectionTool_OpT0FinderApp",
            "-lSelectionTool_OpT0PhotonLibrary","-lSelectionTool_OpT0FinderAlgorithms","-lSelectionTool_OpT0FinderBase"]
    lllibs = [x.replace("lib","-l").replace(".so","") for x in os.listdir(os.environ["LARLITE_LIBDIR"]) if x.endswith('.so') and "LArOpenCV" in x ]
    #libs += ["-L%s"%(os.environ["LARLITE_LIBDIR"])]
    libs += lllibs
if 'GEO2D_BASEDIR' in os.environ:
    libs+= commands.getoutput('geo2d-config --libs').split()
if 'LARCV_BASEDIR' in os.environ:
    libs+= commands.getoutput('larcv-config --libs').split()
if 'ANN_LIBDIR' in os.environ:
    libs+= ["-L%s -lANN" % ( os.environ["ANN_LIBDIR"].strip() )]
    
objs_list=[]
dict_list=[]
for l in dirs:
    d='%s/%s' % (os.environ['LARLITECV_BUILDDIR'],l)
    #print d, os.path.isdir(d)
    #print os.listdir(d)
    src_obj = ['%s/%s' % (d,x) for x in os.listdir(d) if x.endswith('.o') and not x.endswith('Dict.o')]
    dic_obj = ['%s/%s' % (d,x) for x in os.listdir(d) if x.endswith('Dict.o')]

    if len(src_obj) == 0 and len(dic_obj) == 0: continue
    
    objs_list.append(src_obj)
    dict_list.append(dic_obj)
    #print objs_list[-1]
    #print dict_list[-1]

if not dict_list and not objs_list: sys.exit(0)

libname = '%s/liblarlitecv.so' % os.environ['LARLITECV_LIBDIR']
    
cmd='-o %s ' % libname
base_ts = None
if os.path.isfile(libname):
    base_ts=os.path.getmtime(libname)

skip_build = base_ts is not None
for objs in objs_list:

    for obj in objs:

        cmd += '%s ' % obj

        ts = os.path.getmtime(obj)

        if skip_build and ts > base_ts:
            skip_build = False

if skip_build: sys.exit(0)
    
for d in dict_list:

    if d: cmd += '%s ' % d[0]

for l in libs:
    cmd += '%s ' % l

if 'build' in sys.argv: print 1
else: print cmd
sys.exit(1)
