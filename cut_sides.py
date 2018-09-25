import sys
import os
import ntpath


def check_files(file_list):
    with open(file_list) as f:
        for l in f:
            l = l.strip()
            if not os.path.exists(l):
                print "In", (file_list + ","), l, "does not exists..."
                return False
    return True


def cut_side(k, fname, kmc_dir):
    dir_path, fn = ntpath.split(fname)
    begin_fname = os.path.join(kmc_dir, fn + "_k" + str(k - 1) + ".begin")
    end_fname = os.path.join(kmc_dir, fn + "_k" + str(k - 1) + ".end")

    with open(begin_fname, "w") as begin:
        with open(end_fname, "w") as end:
            with open(fname, "r") as f:
                for l in f:
                    if l[0] == ">" or l[1] == ">":
                        begin.write(l)
                        end.write(l)
                    else:
                        begin.write(l[:(k - 1)])
                        begin.write("\n")
                        
                        if l[-1] == "\n":
                            end.write(l[-k:])
                        else:
                            end.write(l[-k + 1 :])


def cut_sides(k, file_list, kmc_dir):
    if k < 2:
        print "k must be >= 2..."
    if not os.path.isdir(kmc_dir):
        print kmc_dir, "does not exists..."
    elif not os.path.exists(file_list):
        print file_list, "does not exists..."
    elif check_files(file_list):
        with open(file_list) as f:
            for l in f:
                l = l.strip()
                cut_side(k, l, kmc_dir)


if __name__ == "__main__":
    if "-k=" not in sys.argv[1] or len(sys.argv) < 4:
        print "Bad argument list..."
        print "example usage: ./colorgra_tool -k=32 file_list.txt kmc-files-dir"
    else:
        cut_sides(int(sys.argv[1][3:]), sys.argv[2], sys.argv[3])
