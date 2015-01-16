//
//  fop.cpp
//  readvtklis
//
//  Created by Rixin Li on 1/15/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include "fop.h"

FileIO::FileIO(int argc, const char * argv[])
{
    int temp;
    int iflag = 0, bflag = 0, sflag = 0, fflag = 0, oflag = 0;
    if (argc < 11) {
        cout << "USAGE: " << argv[0] << " -i <data_path> -b <data_basename> -s <post_name> -f <# (range(f1:f2))> -o <output_path_name>\n" << endl;
        cout << "Example: ./readvtklis -i comb -b Cout -s all -f 0:100 -o result.txt" << endl;
        exit(1);
    } else {
        while ((temp = getopt(argc, (char **)argv, "i:b:s:f:o:")) != -1) {
            switch (temp) {
                case 'i':
                    data_path.assign(optarg);
                    cout << "data_path is " << data_path << endl;
                    iflag = 1;
                    break;
                case 'b':
                    data_basename.assign(optarg);
                    cout << "data_basename is " << data_basename << endl;
                    bflag = 1;
                    break;
                case 's':
                    post_name.assign(optarg);
                    cout << "post_name is " << post_name << endl;
                    sflag = 1;
                    break;
                case 'f':
                    sscanf(optarg,"%d:%d", &start_no, &end_no);
                    if (start_no < 0) {
                        cout << "The start number should be larger than 0. (Auto fix to 0)" << endl;
                        start_no = 0;
                    }
                    if (end_no < start_no) {
                        cout << "The end number should be larger than the start number. (Auto fix to start number + 1)." << endl;
                        end_no += start_no;
                    }
                    cout << "start_no=" << start_no << ", end_no=" << end_no << endl;
                    fflag = 1;
                    break;
                case 'o':
                    output_path_name.assign(optarg);
                    cout << "output_path_name is " << output_path_name;
                    oflag = 1;
                    break;
                case '?':
                    if (optopt == 'i' || optopt == 'b' || optopt == 's' || optopt == 'f' || optopt == 'o')
                        fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                    else if (isprint (optopt))
                        fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                    else
                        fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                    exit(2);
                default:
                    cout << "Argument wrong." << endl;
                    abort();
            }
        }
        // Since argc needs to > 11, so the five IF statements should not be called at all.
        if (iflag == 0) {
            cout << "We need data_path. Abort." << endl;
            abort();
        }
        if (bflag == 0) {
            cout << "We need data_basename. Abort." << endl;
            abort();
        }
        if (sflag == 0) {
            cout << "We need post_name. Abort." << endl;
            abort();
        }
        if (fflag == 0) {
            cout << "We need file range. Abort." << endl;
            abort();
        }
        if (oflag == 0) {
            cout << "We need output_path_name. Abort." << endl;
            abort();
        }
    }
}

int FileIO::Generate_Filename()
{
    
    return 0;
}

