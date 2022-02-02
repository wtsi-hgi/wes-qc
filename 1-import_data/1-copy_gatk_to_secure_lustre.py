#copy GATK variant calls from cromwell execution directory to secure lustre

import os
import subprocess

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
    shard_intervals = {}
    split_interval_dir = wdir + "/call-SplitIntervalList/"
    exec_dir = find_exec_dir(split_interval_dir)
    scatter_dir = exec_dir + "scatterDir/"
    
    scatter_dir_list = os.listdir(scatter_dir)

    for intervalfile in scatter_dir_list:
        intervalfilepath = scatter_dir_list + intervalfile
        intervalnum = intervalfile[0:4]
        #remove leading zeros from file name to get shard number, except for 0
        if intervalnum == '0000':
            intervalnum = '0'
        else:
            while intervalnum[0] == '0':
                intervalnum = intervalnum[1:]

        first_interval_cmd = "grep -m1 ^chr " + intervalfile
        first_interval_output = subprocess.getoutput(first_interval_cmd)
        startchrom = first_interval_output.split()[0]
        startpos = first_interval_output.split()[1]

        last_interval_cmd = "tail -n 1 " + intervalfile
        last_interval_output = subprocess.getoutput(last_interval_cmd)
        endchrom = last_interval_output.split()[0]#start and end of intervals could be on different chrom
        endpos = last_interval_output.split()[2]

        shard_intervals[intervalnum] = ('_').join([startchrom, startpos, endchrom, endpos])

    return shard_intervals

def copy_output_vcfs(wdir, shard_intervals, targetdir):
    '''
    for each shard interval, find the output vcf.gz and tbi in call-MergeVCF and copy to target dir, renaming with intervals
    '''
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
        cp_vcf_cmd = "cp " + vcf_path + " " + dest_vcf
        cp_tabix_cmd = "cp " + tabix_path + " " + dest_tabix
        print(cp_vcf_cmd)
        print(cp_tabix_cmd)
        exit(0)


def main():
    cromwell_base_dir = '/lustre/scratch119/realdata/mdt3/projects/birth_cohort_wes/alspac/gatk/joint_call/cromwell/cromwell-executions/JointCalling/'
    target_dir = '/lustre/scratch123/hgi/mdt2/projects/birth_cohort_wes/qc/gatk_vcfs/'
    regions = {'chr1_3_dir': '303af8c0-4b7b-4388-98b5-2f552584a8f7',
               'chr4_8_dir': '5b20ad2f-9f35-4713-b171-a7f8f0b0aebc',
               'chr9_14_dir': '679738e5-0a41-4196-8b4e-f1e01d80b077',
               'chr15_Y_dir': 'e5540273-472e-4a46-b691-568316b4075d'
               }

    for reg in regions.keys():
        wdir = cromwell_base_dir + regions[reg]
        shard_intervals = get_shard_intervals(wdir)
        copy_output_vcfs(wdir, shard_intervals, target_dir)

if __name__ == '__main__':
    main() 
