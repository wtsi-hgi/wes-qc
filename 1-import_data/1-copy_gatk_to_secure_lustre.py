#copy GATK variant calls from cromwell execution directory to secure lustre

import os
import subprocess
import shutil
import sys
import yaml

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def find_exec_dir(indir):
    '''
    Find most recent execution dir
    '''
    dirlist = os.listdir(indir)
    
    #look for repeated attempts
    repeats = [x for x in dirlist if x.startswith('attempt-')]
    if len(repeats) > 0:#there are repeated dirs
        repnum = len(repeats)
        #if there is one more attempt it is attempt-2, 2 more = attempt-3
        #we want the most recent
        repnum += 1
        exec_dir = indir + "attempt-" + str(repnum) + "/execution/"
        ## TODO FIX THIS TO READ CACHED COPIES ##
#    elif 'cacheCopy' in dirlist:
#        #this assumes cacheDir has no subdirectories for repeat attempts
#        exec_dir = indir + "cacheCopy/execution/"
    elif 'execution' in dirlist:
        exec_dir = indir + "execution/"
    else:
        print("No execution dir found for " + indir)
        return None
    
    if os.path.isdir(exec_dir):
        return exec_dir
    else:
        print("Not a directory " + exec_dir)
        return None


def get_shard_intervals(wdir):
    '''
    create a dict containing each shard number and it's start and end locations
    '''
    print("Creating intervals dict")
    shard_intervals = {}
    split_interval_dir = wdir + "/call-SplitIntervalList/"
    exec_dir = find_exec_dir(split_interval_dir)
    scatter_dir = exec_dir + "scatterDir/"
    
    scatter_dir_list = os.listdir(scatter_dir)

    for intervalfile in scatter_dir_list:
        intervalfilepath = scatter_dir + intervalfile
        intervalnum = intervalfile[0:4]
        #remove leading zeros from file name to get shard number, except for 0
        if intervalnum == '0000':
            intervalnum = '0'
        else:
            while intervalnum[0] == '0':
                intervalnum = intervalnum[1:]

        first_interval_cmd = "grep -m1 ^chr " + intervalfilepath
        first_interval_output = subprocess.getoutput(first_interval_cmd)
        startchrom = first_interval_output.split()[0]
        startpos = first_interval_output.split()[1]

        last_interval_cmd = "tail -n 1 " + intervalfilepath
        last_interval_output = subprocess.getoutput(last_interval_cmd)
        endchrom = last_interval_output.split()[0]#start and end of intervals could be on different chrom
        endpos = last_interval_output.split()[2]

        shard_intervals[intervalnum] = ('_').join([startchrom, startpos, endchrom, endpos])

    return shard_intervals


def copy_output_vcfs(wdir, shard_intervals, targetdir):
    '''
    for each shard interval, find the output vcf.gz and tbi in call-MergeVCF and copy to target dir, renaming with intervals
    '''
    shardcount = len(shard_intervals.keys())
    print("copying and renaming " + str(shardcount) + " VCFs")
    merge_dir = wdir + "/call-MergeVCFs/"

    for shard in shard_intervals.keys():
        shard_dir = merge_dir + "shard-" + shard + "/"
        exec_dir = find_exec_dir(shard_dir)
        vcf_path = exec_dir + 'merged_output.vcf.gz'
        tabix_path = exec_dir + 'merged_output.vcf.gz.tbi'
        #check tabix file exists as this is produced last - quit if it doesn't
        if not os.path.isfile(tabix_path):
            print("Output lacking for " + exec_dir)
            exit(1)
        dest_vcf = targetdir + shard_intervals[shard] + ".gatk.vcf.gz"
        dest_tabix = targetdir +  shard_intervals[shard] + ".gatk.vcf.gz.tbi"
        shutil.copyfile(vcf_path, dest_vcf)
        shutil.copyfile(tabix_path, dest_tabix)


def main():

    script_dir = get_script_path()
    input_yaml = script_dir + '/../config/inputs.yaml'
    with open(input_yaml, 'r') as y:
        inputs = yaml.load(y, Loader=yaml.FullLoader)

    cromwell_base_dir = inputs['cromwell_base_dir']
    target_dir = inputs['target_dir']
    regions = inputs['regions']

    for reg in regions.keys():
        print("Processing " + reg)
        wdir = cromwell_base_dir + regions[reg]
        shard_intervals = get_shard_intervals(wdir)
        copy_output_vcfs(wdir, shard_intervals, target_dir)


if __name__ == '__main__':
    main() 
