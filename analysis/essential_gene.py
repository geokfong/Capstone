import csv

input_file = "CRISPRGeneEffect.csv"
transposed_file = "transposed_CRISPRGeneEffect.csv"
essential_genes_file = "essential_genes.txt"

with open(input_file, newline='') as infile:
    reader = list(csv.reader(infile))
    transposed = list(zip(*reader))  # Transpose

with open(transposed_file, "w", newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerows(transposed)

with open(transposed_file, newline='') as infile, open(essential_genes_file, "w") as outfile:
    reader = csv.reader(infile)
    header = next(reader)

    for row in reader:
        gene = row[0]
        values = []
        for val in row[1:]:
            try:
                values.append(float(val))
            except ValueError:
                continue

        if values:
            avg = sum(values) / len(values)
            if avg < -0.5:
                outfile.write(gene + "\n")
