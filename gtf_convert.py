import re
import csv

def parse_gtf_file(file_path):
    gene_data = {}

    with open(file_path, 'r') as gtf_file:
        for line in gtf_file:
            if line.startswith('#'):
                continue

            attributes = re.findall(r'(\w+)\s+\"(.*?)\"', line)
            record = dict(attributes)

            gene_id = record.get('gene_id')
            gene_name = record.get('gene_name')

            if gene_id and gene_name:
                if gene_id in gene_data:
                    continue

                gene_data[gene_id] = {'gene_id': gene_id, 'gene_name': gene_name}

    return list(gene_data.values())

def write_to_csv(data, output_file):
    with open(output_file, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['Gene ID', 'Gene Name'])
        for record in data:
            writer.writerow([record['gene_id'], record['gene_name']])

parsed_data = parse_gtf_file("gencode.vM32.primary_assembly.annotation.gtf")

print(parsed_data)
write_to_csv(parsed_data, 'gene_data.csv')
