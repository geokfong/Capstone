echo "Creating comprehensive Excel-ready output file..."
awk -F'\t' 'NR>1 {print $1 "\t" $2 "\t" $3 "\t" $4}' bt_sample_edits/feature_region_edit_spans.bed > bt_sample_edits/temp_coordinates.txt
{
    echo -e "chrom\tstart\tend\tcluster_info\tscore\tstrand\tregion_length\tgwas_diseases\tgwas_count\tclinvar_diseases\tclinvar_count\ttotal_disease_count\tdisease_sources\tfeature_type\trepeat_elements\trepeat_count\trepeat_types\trepeat_coverage_bp\thighest_expression_tissue\thighest_expression_value\tpriority_score\tgenomic_location"
    while IFS=$'\t' read -r chrom start end cluster_info; do
        region_line=$(awk -F'\t' -v c="$chrom" -v s="$start" -v e="$end" -v ci="$cluster_info" \
            '$1==c && $2==s && $3==e && $4==ci {print; exit}' bt_sample_edits/feature_region_edit_spans.bed)
        
        if [ -n "$region_line" ]; then
            score=$(echo "$region_line" | cut -f5)
            strand=$(echo "$region_line" | cut -f6)
            region_length=$((end - start))
            genomic_location="${chrom}:${start}-${end}"
            
            feature_type="other"
            feature_priority=5
            case "$cluster_info" in
                *intron*) feature_type="intron"; feature_priority=1 ;;
                *3_prime_UTR*|*3UTR*) feature_type="3_UTR"; feature_priority=2 ;;
                *5_prime_UTR*|*5UTR*) feature_type="5_UTR"; feature_priority=3 ;;
                *exon*) feature_type="exon"; feature_priority=4 ;;
            esac
            
            # Get repeat elements from RepeatMasker data 
            repeat_info=$(awk -F'\t' -v c="$chrom" -v s="$start" -v e="$end" -v ci="$cluster_info" \
                'BEGIN {
                    repeat_types = ""
                    repeat_count = 0
                    total_coverage = 0
                    delete seen_types  # Clear array
                }
                $1==c && $2==s && $3==e && $4==ci && $10!="." && $10!="" {
                    # This region has repeat overlaps
                    repeat_type = $10
                    repeat_length = ($12 > 0) ? $12 : ($9 - $8)  # Use overlap length or repeat length
                    
                    # Track unique repeat types
                    if (seen_types[repeat_type] == 0) {
                        seen_types[repeat_type] = 1
                        if (repeat_types == "") {
                            repeat_types = repeat_type
                        } else {
                            repeat_types = repeat_types "," repeat_type
                        }
                    }
                    repeat_count++
                    total_coverage += repeat_length
                }
                END {
                    if (repeat_types == "") {
                        repeat_types = "none"
                        repeat_count = 0
                        total_coverage = 0
                        elements = "none"
                    } else {
                        elements = "present"
                    }
                    print elements "\t" repeat_count "\t" repeat_types "\t" total_coverage
                }' bt_sample_edits/output_repeatmasker_class.bed)
            
            repeat_elements=$(echo "$repeat_info" | cut -f1)
            repeat_count=$(echo "$repeat_info" | cut -f2)
            repeat_types=$(echo "$repeat_info" | cut -f3)
            repeat_coverage=$(echo "$repeat_info" | cut -f4)
            
            [ -z "$repeat_elements" ] && repeat_elements="none"
            [ -z "$repeat_count" ] && repeat_count=0
            [ -z "$repeat_types" ] && repeat_types="none"
            [ -z "$repeat_coverage" ] && repeat_coverage=0
            
            # Get GWAS diseases 
            gwas_info=$(awk -F'\t' -v c="$chrom" -v s="$start" -v e="$end" -v ci="$cluster_info" \
                'BEGIN {diseases=""; count=0}
                $1==c && $2==s && $3==e && $4==ci && $16!="" && $16!="." && $16!="NA" {
                    # Clean and sanitize disease name
                    disease = $16
                    gsub(/[\t\n\r]/, " ", disease)  # Replace tabs/newlines with spaces
                    gsub(/;/, ",", disease)         # Replace semicolons with commas
                    gsub(/  +/, " ", disease)       # Replace multiple spaces with single space
                    gsub(/^ +| +$/, "", disease)    # Trim leading/trailing spaces
                    
                    if (seen[disease] == 0) {
                        seen[disease] = 1
                        if (diseases == "") {
                            diseases = disease
                        } else {
                            diseases = diseases ";" disease
                        }
                        count++
                    }
                } 
                END {
                    if (diseases == "") diseases = "none"
                    print diseases "\t" count
                }' bt_sample_edits/output_GWAS.tsv)
                
            gwas_diseases=$(echo "$gwas_info" | cut -f1)
            gwas_count=$(echo "$gwas_info" | cut -f2)
            [ -z "$gwas_diseases" ] && gwas_diseases="none" && gwas_count=0
            
            # Get ClinVar diseases
            clinvar_info=$(awk -F'\t' -v c="$chrom" -v s="$start" -v e="$end" -v ci="$cluster_info" \
                'BEGIN {diseases=""; count=0}
                $1==c && $2==s && $3==e && $4==ci && $29!="" && $29!="." && $29!="NA" {
                    # Clean and sanitize disease name
                    disease = $29
                    gsub(/[\t\n\r]/, " ", disease)  # Replace tabs/newlines with spaces
                    gsub(/;/, ",", disease)         # Replace semicolons with commas
                    gsub(/  +/, " ", disease)       # Replace multiple spaces with single space
                    gsub(/^ +| +$/, "", disease)    # Trim leading/trailing spaces
                    
                    if (seen[disease] == 0) {
                        seen[disease] = 1
                        if (diseases == "") {
                            diseases = disease
                        } else {
                            diseases = diseases ";" disease
                        }
                        count++
                    }
                } 
                END {
                    if (diseases == "") diseases = "none"
                    print diseases "\t" count
                }' bt_sample_edits/output_clinvarMain.tsv)
                
            clinvar_diseases=$(echo "$clinvar_info" | cut -f1)
            clinvar_count=$(echo "$clinvar_info" | cut -f2)
            [ -z "$clinvar_diseases" ] && clinvar_diseases="none" && clinvar_count=0
            
            # Calculate disease info
            total_disease_count=$((gwas_count + clinvar_count))
            
            # Determine disease sources and priority
            if [ "$gwas_count" -gt 0 ] && [ "$clinvar_count" -gt 0 ]; then
                disease_sources="Both"
                source_priority=1
            elif [ "$clinvar_count" -gt 0 ]; then
                disease_sources="ClinVar"
                source_priority=2
            elif [ "$gwas_count" -gt 0 ]; then
                disease_sources="GWAS"
                source_priority=3
            else
                disease_sources="none"
                source_priority=4
            fi
            
            priority_score=$((source_priority * 10 + feature_priority))
            
            highest_expression_tissue="no_data"
            highest_expression_value=0
            
            if [ -f "bt_sample_edits/tissue_expression/combined_expression_matrix.tsv" ]; then
                expr_info=$(awk -F'\t' -v c="$chrom" -v s="$start" -v e="$end" \
                    'NR==1 {
                        # Store tissue names from header
                        for (i=8; i<=NF; i+=3) {
                            tissue_idx = (i-5)/3
                            tissue_name[tissue_idx] = $i
                            gsub(/_mean$/, "", tissue_name[tissue_idx])
                        }
                    }
                    NR>1 && $1==c && $2==s && $3==e {
                        max_expr = 0
                        max_tissue = "no_expression"
                        
                        for (i=8; i<=NF; i+=3) {
                            if ($i != "" && $i != "." && $i > max_expr) {
                                max_expr = $i
                                tissue_idx = (i-5)/3
                                max_tissue = tissue_name[tissue_idx]
                            }
                        }
                        
                        if (max_tissue == "") max_tissue = "no_expression"
                        print max_tissue "\t" max_expr
                        exit
                    }' bt_sample_edits/tissue_expression/combined_expression_matrix.tsv)
                
                if [ -n "$expr_info" ]; then
                    highest_expression_tissue=$(echo "$expr_info" | cut -f1)
                    highest_expression_value=$(echo "$expr_info" | cut -f2)
                fi
            fi
            
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$chrom" "$start" "$end" "$cluster_info" "$score" "$strand" "$region_length" \
                "$gwas_diseases" "$gwas_count" "$clinvar_diseases" "$clinvar_count" \
                "$total_disease_count" "$disease_sources" "$feature_type" "$repeat_elements" \
                "$repeat_count" "$repeat_types" "$repeat_coverage" "$highest_expression_tissue" \
                "$highest_expression_value" "$priority_score" "$genomic_location"
        fi
    done < bt_sample_edits/temp_coordinates.txt
    
} | sort -t$'\t' -k21,21n -k1,1 -k2,2n > bt_sample_edits/comprehensive_analysis_results_excel.tsv
# Clean up temporary file
rm -f bt_sample_edits/temp_coordinates.txt
