write_xml <- function(output_file, 
                      dir_results, mat_proba_migration,
                      sequence_length = 3000,
                      mutation_rate = 0.00005,
                      beta_rate = 0.0001,
                      rate_out_of_E = 0.33,
                      rate_out_of_I = 0.33,
                      p_sample = 0.1,
                      n_demes = 3,
                      min_sample_per_deme = 50,
                      vec_S_init_per_deme = c(1e5, 1e5, 1e5)){
  cat('<beast version="2.0" namespace="beast.base.inference.parameter:beast.base.inference:remaster">','\n',
      file = output_file, sep = '')
  cat('\t', file = output_file, sep = '', append = T)
  cat('<run spec="Simulator" nSims="1">', '\n', 
      file = output_file, sep = '', append = T)
  cat('\t\t', file = output_file, sep = '', append = T)
  cat('<simulate spec="feast.simulation.SimulatedAlignment" ', 
      'outputFileName="', dir_results,
      '/$(filebase).nexus" sequenceLength="', format(sequence_length, scientific = F), 
      '">', '\n',
      file = output_file, sep = '', append = T)
  cat('\t\t\t', file = output_file, sep = '', append = T)
  cat('<siteModel spec="beast.base.evolution.sitemodel.SiteModel" mutationRate="', format(mutation_rate, scientific = F), '">', '\n',
      file = output_file, sep = '', append = T)
  cat('\t\t\t\t', file = output_file, sep = '', append = T)
  cat('<substModel spec="beast.base.evolution.substitutionmodel.JukesCantor"/>', '\n', 
      file = output_file, sep = '', append = T)
  cat('\t\t\t', file = output_file, sep = '', append = T)
  cat('</siteModel>', '\n', 
      file = output_file, sep = '', append = T)
  cat('\n\n', file = output_file, sep = '', append = T)
  cat('\t\t\t', file = output_file, sep = '', append = T)
  
  cat('<tree spec="SimulatedTree" id="tree">', '\n', file = output_file, sep = '', append = T)
  cat('\t\t\t\t', file = output_file, sep = '', append = T)
  cat('<trajectory spec="StochasticTrajectory" id="traj" mustHave="', file = output_file, sep = '', append = T)
  for(deme in 0:(n_demes - 1)){
    cat('sample[', deme, ']&gt;', min_sample_per_deme, sep = '', file = output_file, append = T)
    if(deme != n_demes - 1){
      cat(' &amp;&amp; ', sep = '', file = output_file, append = T)
    }
  }
  cat('">', '\n\n', file = output_file, sep = '', append = T)
  
  ## Population initialization
  cat('\t\t\t\t\t', file = output_file, sep = '', append = T)
  cat('<population spec="RealParameter" id="S" value="', file = output_file, sep = '', append = T)
  cat(format(vec_S_init_per_deme, scientific = F), file = output_file,append = T)
  cat('"/>', '\n', file = output_file, sep = '', append = T)
  
  cat('\t\t\t\t\t', file = output_file, sep = '', append = T)
  cat('<population spec="RealParameter" id="E" value="', file = output_file, sep = '', append = T)
  cat(rep(0, n_demes), file = output_file,append = T)
  cat('"/>', '\n', file = output_file, sep = '', append = T)
  
  cat('\t\t\t\t\t', file = output_file, sep = '', append = T)
  cat('<population spec="RealParameter" id="I" value="', file = output_file, sep = '', append = T)
  cat(c(1, rep(0, n_demes - 1)), file = output_file,append = T)
  cat('"/>', '\n', file = output_file, sep = '', append = T)
  
  cat('\t\t\t\t\t', file = output_file, sep = '', append = T)
  cat('<population spec="RealParameter" id="R" value="', file = output_file, sep = '', append = T)
  cat(rep(0, n_demes), file = output_file,append = T)
  cat('"/>', '\n', file = output_file, sep = '', append = T)
  
  cat('\t\t\t\t\t', file = output_file, sep = '', append = T)
  cat('<samplePopulation spec="RealParameter" id="sample" value="', file = output_file, sep = '', append = T)
  cat(rep(0, n_demes), file = output_file,append = T)
  cat('"/>', '\n', file = output_file, sep = '', append = T)
  cat('\n', file = output_file, sep = '', append = T)
  
  # S to E transition
  for(deme_origin in 0:(n_demes - 1)){
    for(deme_destination in (0:(n_demes - 1))){
      
      cat('\t\t\t\t\t', file = output_file, sep = '', append = T)
      cat('<reaction spec="Reaction" rate="', 
          format(beta_rate * mat_proba_migration[deme_origin + 1, deme_destination + 1], scientific = F),
          '"> ',
          'I[', deme_origin, '] + S[', deme_destination, '] -> I[', deme_origin, '] + E[', deme_destination, '] </reaction>',
          '\n',
          file = output_file, sep = '', append = T)
    }
  }
  cat('\n', file = output_file, sep = '', append = T)
  
  # E to I transitions
  for(deme in 0:(n_demes - 1)){
    cat('\t\t\t\t\t', file = output_file, sep = '', append = T)
    cat('<reaction spec="Reaction" rate="', format(rate_out_of_E, scientific = F), '"> ',
        'E[', deme, '] -> I[', deme, '] </reaction>',
        '\n',
        file = output_file, sep = '', append = T)
  }
  cat('\n', file = output_file, sep = '', append = T)
  
  # I to R transitions
  for(deme in 0:(n_demes - 1)){
    cat('\t\t\t\t\t', file = output_file, sep = '', append = T)
    cat('<reaction spec="Reaction" rate="', format(rate_out_of_I * (1. - p_sample), scientific = F), '"> ',
        'I[', deme, '] -> R[', deme, '] </reaction>',
        '\n',
        file = output_file, sep = '', append = T)
  }
  cat('\n', file = output_file, sep = '', append = T)
  
  # I to sample transitions
  for(deme in 0:(n_demes - 1)){
    cat('\t\t\t\t\t', file = output_file, sep = '', append = T)
    cat('<reaction spec="Reaction" rate="', format(rate_out_of_I * p_sample, scientific = F), '"> ',
        'I[', deme, '] -> sample[', deme, '] </reaction>',
        '\n',
        file = output_file, sep = '', append = T)
  }
  cat('\n', file = output_file, sep = '', append = T)
  
  cat('\t\t\t\t', '</trajectory>', '\n', file = output_file, sep = '', append = T)
  cat('\t\t\t', '</tree>', '\n', file = output_file, sep = '', append = T)
  cat('\t\t', '</simulate>', '\n\n', file = output_file, sep = '', append = T)
  
  cat('\t\t', '<logger spec="Logger" fileName="',
      dir_results, 
      '/$(filebase).traj">', '\n', file = output_file, sep = '', append = T)
  cat('\t\t\t', '<log idref="traj"/>', '\n', file = output_file, sep = '', append = T)
  cat('\t\t', '</logger>', '\n\n', file = output_file, sep = '', append = T)
  
  cat('\t\t', ' <logger spec="Logger" mode="tree" fileName="',
      dir_results, '/$(filebase).trees">', '\n', file = output_file, sep = '', append = T)
  cat('\t\t\t', '<log spec="TypedTreeLogger" tree="@tree" removeSingletonNodes="true"/>', '\n', file = output_file, sep = '', append = T)
  cat('\t\t', '</logger>', '\n\n', file = output_file, sep = '', append = T)
  
  cat('\t\t', ' <logger spec="Logger">', '\n', file = output_file, sep = '', append = T)
  cat('\t\t\t', '<log spec="beast.base.evolution.tree.TreeStatLogger" tree="@tree"/>', '\n', file = output_file, sep = '', append = T)
  cat('\t\t', '</logger>', '\n\n', file = output_file, sep = '', append = T)
  
  cat('\t', '</run>', '\n', file = output_file, sep = '', append = T)
  cat('</beast>', '\n', file = output_file, sep = '', append = T)
}

simulate_from_xml <- function(beast_path, xml_file_path){
  system(paste0(beast_path, ' -overwrite ', xml_file_path))
}

parse_nexus_remaster_save_fasta <- function(nexus_file_path, fasta_file_path){
  con <- file(nexus_file_path, "r")
  
  pre_sequences <- readLines(con, n = 11)
  n_sequences <- strsplit(strsplit(x = pre_sequences[4], split = '=')[[1]][2], split = ';')[[1]][1] %>% 
    as.numeric()
  
  list_name_seq <- vector('list', n_sequences)
  list_seq <- vector('list', n_sequences)
  
  for(i_seq in 1:n_sequences){
    line_seq <- readLines(con, n = 1)
    line_seq <- substr(x = line_seq, start = 3, stop = nchar(line_seq))
    line_seq <- strsplit(line_seq, ' ')[[1]]
    list_name_seq[[i_seq]] <- line_seq[1]
    list_seq[[i_seq]] <- strsplit(line_seq[2], split = ';')[[1]][1]
  }
  
  close(con)
  
  names(list_seq) <- Reduce('c', list_name_seq)
  
  write.FASTA(char2dna(list_seq), file = fasta_file_path)
}

# Function to save metadata file from the phylo_tree output from remaster
write_metadata_from_phylo_tree <- function(metadata_file_path, phylo_tree){
  metadata_seq <- phylo_tree %>% 
    filter(isTip == T, !is.na(time)) %>%
    rename(subgroup = type) %>% 
    mutate(subgroup = substr(subgroup, start = 3, stop = nchar(subgroup))) %>% 
    select(label, time, subgroup) %>% 
    as_tibble() %>% 
    filter(! is.na(subgroup)) %>% 
    rename(strain = label)
  
  write.csv(metadata_seq, metadata_file_path, quote = F, row.names = F)
}

get_alignment_metadata <- function(xml_file_path, dir_results, name_scenario){
  # Read transmission tree (with associated node characteristics)
  phylo_tree <- read.beast(file = paste0(dir_results, '/', scenario_name, '.trees'))
  
  # Save alignment data from nexus
  parse_nexus_remaster_save_fasta(paste0(dir_results, '/', scenario_name, '.nexus'), 
                                  paste0(dir_results, '/alignment_', name_scenario, '.fasta'))
  
  write_metadata_from_phylo_tree(paste0(dir_results, '/metadata_', name_scenario, '.csv'), phylo_tree)

}
