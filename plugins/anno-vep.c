/*  plugins/anno-vep.c -- adds tags to the CSQ field in VEP-annotated vcfs

    Copyright (C) 2022 Marco Baggio

    Author: Marco Baggio <marco.baggio@kuleuven.be>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <htslib/vcf.h>
#include "../cols.h"

typedef struct {
    char* key;
    char* value;
} item;

int cmp(const void* a, const void* b) {
    item* item_a = (item*)a;
    item* item_b = (item*)b;
    return strcmp(item_a->key, item_b->key);
}

int parse_file(const char* filename, item** il_pt) {
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) return 0;

    int alloc_steps = 1000;
    int alloc = alloc_steps;
    int nrecords = 0;
    *il_pt = malloc(alloc * sizeof(item));

    char* line = malloc(1024 * sizeof(char));
    size_t n;
    int line_len = 0;
    char* key;
    char* value;

    while (getline(&line, &n, fp) > 0) {
        line_len = strlen(line);
        if (line[line_len - 1] == '\n') line[line_len - 1] = 0;
        key = strtok(line, "\t");
        value = strtok(NULL, "\t");

        if (key && value) {
            nrecords++;
            if (nrecords > alloc) {
                alloc += alloc_steps;
                *il_pt = realloc(*il_pt, alloc * sizeof(item));
            }
            (*il_pt)[nrecords - 1].key = strdup(key);
            (*il_pt)[nrecords - 1].value = strdup(value);
        }
    }

    fclose(fp);
    return nrecords;
}

bcf_hdr_t* header;
char* csq_str;
int ncsq_str;
cols_t* cols_tr,    // the current CSQ tag split into transcripts
* cols_csq;         // the current CSQ transcript split into fields

item* items;
int nitems;

/*
    This short description is used to generate the output of `bcftools plugin -l`.
*/
const char* about(void) {
    return
        "A plugin to add tags to the CSQ field in VEP-annotated VCFs\n"
        "Usage: bcftools +anno-vep <in.vcf> -- TAG_NAME TSV_FILE\n";
}

/*
    Called once at startup, allows to initialize local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/
int init(int argc, char** argv, bcf_hdr_t* in, bcf_hdr_t* out) {
    if (argc < 3) {
        printf("%s", about());
        return -1;
    }

    char* new_tag = strdup(argv[1]);
    char* filename = strdup(argv[2]);

    nitems = parse_file(filename, &items);
    if (nitems < 1) {
        printf("Error reading the file %s\n",filename);
        return -1;
    }

    header = out;
    bcf_hrec_t* hrec = bcf_hdr_get_hrec(header, BCF_HL_INFO, NULL, "CSQ", NULL);
    if (hrec) {
        int ret = bcf_hrec_find_key(hrec, "Description");
        if (ret >= 0) {
            char* format = strstr(hrec->vals[ret], "Format: ");
            if (format) {
                format += 8;
                format[strlen(format) - 1] = 0;
            }
            char* new_ptr = realloc(hrec->vals[ret], strlen(hrec->vals[ret]) + strlen(new_tag) + 3);
            if (new_ptr) {
                hrec->vals[ret] = new_ptr;
                strcat(hrec->vals[ret], "|");
                strcat(hrec->vals[ret], new_tag);
                strcat(hrec->vals[ret],"\"");
            } else {
                printf("Out of memory\n");
                return -1;
            }
        }
    }
    return 0;
}


/*
    Called for each VCF record. Return rec to output the line or NULL
    to suppress output.
*/
bcf1_t* process(bcf1_t* rec) {
    int lcsq_str = bcf_get_info_string(header, rec, "CSQ", &csq_str, &ncsq_str);

    if (lcsq_str > 0) {
        char* gene_id;
        char* transcript;
        item key;
        key.value = strdup("");

        cols_tr = cols_split(csq_str, cols_tr, ',');
        char* new_csq_str = calloc(1, strlen(csq_str) + 1024 * cols_tr->n);

        int i;
        for (i = 0; i < cols_tr->n; i++) {
            transcript = cols_tr->off[i];
            strcat(new_csq_str, transcript);
            strcat(new_csq_str, "|");
            cols_csq = cols_split(transcript, cols_csq, '|');
            gene_id = cols_csq->off[4];
            if (strlen(gene_id) > 0) {
                key.key = strdup(gene_id);
                item* found = bsearch(&key, items, nitems, sizeof(item), cmp);
                if (found)
                    strcat(new_csq_str, found->value);
                free(key.key);
            }
            if (i < cols_tr->n - 1)
                strcat(new_csq_str, ",");
        }
        bcf_update_info_string(header, rec, "CSQ", new_csq_str);
        free(key.value);
        free(new_csq_str);
    }
    return rec;
}


/*
    Clean up.
*/
void destroy(void) {
    free(csq_str);
    free(cols_csq);
    free(cols_tr);
    int i;
    for (i = 0; i < nitems; i++) {
        free(items[i].key);
        free(items[i].value);
    }
    free(items);
}
