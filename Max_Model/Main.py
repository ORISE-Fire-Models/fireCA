import os, glob
import shutil
import subprocess
import sys
import FireCA as fca
import numpy as np

def doCA(indir, outdir, fireID, res, maxIter, n_burnt, plotting, hasBarrier, isSimpleCell, growthWeight):
    print(indir)
    CA = fca.fireCA(indir)
    CA.plotting = plotting
    CA.hasBarrier = hasBarrier
    CA.isSimpleCell = isSimpleCell
    CA.consts['fireID'] = fireID
    CA.consts['gif_dir'] = outdir + '/' + fireID + '/' + str(res) + 'm_' + str(maxIter) + 'hrs_' + str(growthWeight) + 'weight/'
    if isSimpleCell:
        CA.consts['gif_dir'] = CA.consts['gif_dir'] + 'SimpleCell/'
    else:
        CA.consts['gif_dir'] = CA.consts['gif_dir'] + 'NormalCell/'
    CA.consts['n_burnt'] = n_burnt
    CA.consts['growthWeight'] = growthWeight
    if os.path.exists(CA.consts['gif_dir']):
        shutil.rmtree(CA.consts['gif_dir'])
    os.makedirs(CA.consts['gif_dir'])

    CA.consts['spatial_resolution'] = res
    CA.consts['maxTimeStep'] = maxIter

    CA.get_data()
    
    print(f"running {indir}")
    CA.preprocess()
    CA.run()
    CA.saveArrivalMap()
        
    np.savetxt(CA.consts['gif_dir']+'stats.csv', np.vstack((CA.CA_burnt_number, CA.FP, CA.FN, CA.kappa)).T, delimiter=',', header='model, FP, FN, kappa')
        
    print("making gif")
    file_list = []
    ffconcat = ["ffconcat version 1.0\n"]
    os.chdir(CA.consts['gif_dir'])

    for filename in glob.glob('*.png'):
        file_list.append("file " + filename + "\n")
    file_list.sort()
    print(file_list)
    datetime_list = [np.datetime64(str(np.datetime64(int(fname[5:-5].split("_")[0]),'D'))+"T"+fname[5:-5].split("_")[1]) for fname in file_list]
    duration_list = np.append(np.diff(datetime_list).astype('timedelta64[s]'), [3600]).astype(float)*32/3600/240
    duration_list = [f"duration {duration:.3f}\n" for duration in duration_list]
    ffconcat += [item for pair in zip(file_list, duration_list) for item in pair]
    ffconcat += ffconcat[-2:]
    fffilename = 'ffconcat.txt'
    with open(fffilename, 'w') as f:
        f.writelines(ffconcat)
    #subprocess.run(f"ffmpeg -y -f concat -safe 0 -i {fffilename} {CA.consts['spatial_resolution']}m.gif".split(" "))

def handle_args():
    indir=""
    outdir=""
    fireID=""
    res=60
    maxIter=100
    growthWeight=1
    n_burnt=False
    plotting=True
    hasBarrier=False
    isSimpleCell=False
    args = sys.argv[1:]
    for arg in args:
        arg = arg.split("=")
        arg[0] = arg[0].strip()
        arg[0] = arg[0].lower()
        if arg[0] == "--indir":
            indir = arg[1].strip()
        elif arg[0] == "--outdir":
            outdir = arg[1].strip()
        elif arg[0] == "--fireid":
            fireID = arg[1].strip()
        elif arg[0] == "--res":
            res = int(arg[1].strip())
        elif arg[0] == "--max_iteration":
            maxIter = float(arg[1].strip())
        elif arg[0] == "--growthweight":
            growthWeight = float(arg[1].strip())
        elif arg[0] == "--n_burnt":
            if arg[1] == "True":
                n_burnt = True
            elif arg[1] == "False":
                n_burnt = False
            else:
                print("--nburnt must be True or False")
                sys.exit()
        elif arg[0] == "--plotting":
            if arg[1] == "True":
                plotting = True
            elif arg[1] == "False":
                plotting = False
            else:
                print("--plotting must be True or False")
                sys.exit()
        elif arg[0] == "--hasbarrier":
            if arg[1] == "True":
                hasBarrier = True
            elif arg[1] == "False":
                hasBarrier = False
            else:
                print("--hasBarrier must be True or False")
        elif arg[0] == "--issimplecell":
            if arg[1] == "True":
                isSimpleCell = True
            elif arg[1] == "False":
                isSimpleCell = False
            else:
                print("--isSimpleCell must be True or False")
                sys.exit()
        else:
            print("invalid args")
            sys.exit()

    return indir, outdir, fireID, res, maxIter, n_burnt, plotting, hasBarrier, isSimpleCell, growthWeight
    
def main():
    indir, outdir, fireID, res, maxIter, n_burnt, plotting, hasBarrier, isSimpleCell, growthWeight = handle_args()
    doCA(indir, outdir, fireID, res, maxIter, n_burnt, plotting, hasBarrier, isSimpleCell, growthWeight)

if __name__ == "__main__":
    main()
