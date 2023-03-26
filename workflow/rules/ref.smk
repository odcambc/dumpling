rule prepare_dict_GATK:
  input:
    expand("{ref_dir}/{{reference}}.fasta", ref_dir=config["ref_dir"])
  output:
    expand("{ref_dir}/{{reference}}.dict", ref_dir=config["ref_dir"])
  log:
    "logs/gatk/{reference}.create_dict.log"
  shell:
    "gatk CreateSequenceDictionary -R {input}"

rule prepare_index:
  input:
    expand("{ref_dir}/{{reference}}.fasta", ref_dir=config["ref_dir"])
  output:
    expand("{ref_dir}/{{reference}}.fasta.fai", ref_dir=config["ref_dir"])
  log:
    "logs/gatk/{reference}.create_dict.log"
  shell:
    "samtools faidx {input}"
