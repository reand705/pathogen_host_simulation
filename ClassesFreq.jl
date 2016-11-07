#= Classes =#
using StatsBase
# Host --------------------------------------------------------------------------
type Host
  # TODO: Is this V constrain enough?
  sequences::BitArray
  # TODO: Initialize with 0 ?
  number_of_times_that_resisted::Int32
end

function resisted!(host::Host)
  host.number_of_times_that_resisted+=1
end

function reproduce(host1::Host, host2::Host, mutation_rate::Float64)
  # TODO: Improve random haplotype selection
  seq_from_1 = rand(1:2)
  seq_from_2 = rand(1:2)
  return (Host(
  [mutate_haplotype(mutation_rate, host1.sequences[:,seq_from_1]) mutate_haplotype(mutation_rate, host2.sequences[:,seq_from_2])],
  0))
end

# PopHost -----------------------------------------------------------------------
type PopHost
  # TODO: Constrain the structure of the sequence - Static Types
  population_of_hosts::Array{Host, 1}
  mutation_rate::Float64
end

function number_of_alleles(population_of_hosts::PopHost)
  #devolve o numero de alelos distintos da populacao
  length = size(population_of_hosts.population_of_hosts[1].sequences[:,1],1)
  #alleles = falses(2^length)
  frequency = fill(Int32(0), 2^length)
  counters = Int32[]
  for host in population_of_hosts.population_of_hosts
    for i = 1:2
      #alleles[to_Integer(host.sequences[:,i]) + 1] = true
      frequency[to_Integer(host.sequences[:,i]) + 1] += Int32(1)
    end
  end
  #for i = 1:2^length
     #if alleles[i]
       #push!(counters, frequency[i])
     #end
  #end
  #counter = 0
  #for allele in alleles
    #if allele
      #counter+=1
    #end
  #end
  return (frequency)
end

function fitness(population_of_hosts::PopHost)
  #Calculates and return the fitness of a population of hosts
  number_of_hosts = size(population_of_hosts.population_of_hosts,1)
  fit = zeros(number_of_hosts)
  sum_of_squares = 0.0
  #Para cada Host
  for i = 1:number_of_hosts
     #verifique quantas vezes ele resistiu
     fit[i] = 1.0*population_of_hosts.population_of_hosts[i].number_of_times_that_resisted
     fit[i] = fit[i]*fit[i]
     sum_of_squares = sum_of_squares + fit[i]
  end
  #Se a soma for 0, todos tem a mesma fitness (ou prob de ser escolhido)
  if sum_of_squares == 0.0
     for j = 1:number_of_hosts
       fit[j] = 1.0/number_of_hosts
     end
  #do contrario, calcule a fracao
  else
     for j = 1:number_of_hosts
       fit[j] = fit[j]/sum_of_squares
     end
  end
  return fit
end

function random_hosts(population_of_hosts::PopHost, size_of_sample::Int32)
  #Return a sample of "size_of_sample" hosts, picked without replacement and with probability proportional to it's fitness
  fit = fitness(population_of_hosts)
  return(sample(population_of_hosts.population_of_hosts, WeightVec(fit), size_of_sample, replace=false))
end

function reproduce(population_of_hosts::PopHost)
  number_of_hosts = size(population_of_hosts.population_of_hosts,1)
  # TODO: Is it faster to initialize with zeros ou with generic Hosts? How to to do it in the second case?
  new_population = fill(population_of_hosts.population_of_hosts[1], number_of_hosts)
  for i = 1:number_of_hosts
    # TODO: Is it ok to create a new variable here ?
    (host1, host2) = random_hosts(population_of_hosts, Int32(2))
    new_population[i] = reproduce(host1, host2, population_of_hosts.mutation_rate)
  end
  return(PopHost(new_population, population_of_hosts.mutation_rate))
end

# Path --------------------------------------------------------------------------
type Path
  # TODO: Constrain the structure of the sequence - Static Types
  sequences::BitArray
  # TODO: Initialize with 0 ?
  number_of_times_that_infected::Int32
end

function infected!(path::Path)
  path.number_of_times_that_infected+=1
end

function reproduce(path::Path, mutation_rate::Float64)
  return(Path(mutate_haplotype(mutation_rate, path.sequences), 0))
end

function pathogen_infects_host(path::Path, host::Host, minimum_length::Int32)
  number_of_peptides = size(path.sequences, 1)
  for i = 1:number_of_peptides
    if host_presents_peptide(path.sequences[:,i], host, minimum_length)
      return(false)
    end
  end
  return(true)
end

# PopPath -----------------------------------------------------------------------
type PopPath
  # TODO: Constrain the structure of the sequence - Static Types
  population_of_pathogens::Array{Path, 1}
  mutation_rate::Float64
end

function fitness(population_of_pathogens::PopPath)
  #Calculates and return the fitness of a population of pathogens
  number_of_pathogens = size(population_of_pathogens.population_of_pathogens,1)
  fit = zeros(number_of_pathogens)
  sum_of_squares = 0.0
  #Para cada Path
  for i = 1:number_of_pathogens
     #verifique quantos ele infectou
     fit[i] = 1.0*population_of_pathogens.population_of_pathogens[i].number_of_times_that_infected
     fit[i] = fit[i]*fit[i]
     sum_of_squares = sum_of_squares + fit[i]
  end
  #Se a soma for 0, todos tem a mesma fitness
  if sum_of_squares == 0.0
     for j = 1:number_of_pathogens
       fit[j] = 1.0/number_of_pathogens
     end
  #do contrario, calcule a fracao
  else
     for j = 1:number_of_pathogens
       fit[j] = fit[j]/sum_of_squares
     end
  end
  return fit
end

function random_pathogen(population_of_pathogens::PopPath, size_of_sample::Int32)
  #Return a sample of "size_of_sample" hosts, picked without replacement and with probability proportional to it's fitness
  #Returns an array if size_of_sample > 1
  fit = fitness(population_of_pathogens)
  return(sample(population_of_pathogens.population_of_pathogens, WeightVec(fit), size_of_sample, replace=false))
end

function reproduce(population_of_pathogens::PopPath)
  number_of_pathogens = size(population_of_pathogens.population_of_pathogens,1)
  # TODO: Is it faster to initialize with zeros ou with generic Path? How to to do it in the second case?
  new_population = fill(population_of_pathogens.population_of_pathogens[1], number_of_pathogens)
  for i = 1:number_of_pathogens
    # TODO: Is it ok to create a new variable here ?
    path1 = random_pathogen(population_of_pathogens, 1)
    # TODO: Is it ok to just access the first element of the array ?
    new_population[i] = reproduce(path1[1], population_of_pathogens.mutation_rate)
  end
  return(PopPath(new_population, population_of_pathogens.mutation_rate))
end


# Misc --------------------------------------------------------------------------
function mutate_alele(mutation_rate::Float64, alele::BitArray)
  # - Warning - : it transforms the alele!
  # TODO: Constrain the structure of the alele - Static Types
  length = size(alele,1)
  if rand() < mutation_rate
    for j = 1:length
      # TODO: iS THERE A BETTER WAY TO ASSIGN THE NEW ARRAY?
      alele[j] = rand(Bool)
    end
  end
  return(alele)
end

function mutate_haplotype(mutation_rate::Float64, haplotype::BitArray)
  # TODO: Constrain the structure of the haplotype - Static Types
  for i = 1:size(haplotype,2)
    # TODO: Is there a better way to select, alocate and return arrays?
    haplotype[:,i] = mutate_alele(mutation_rate, haplotype[:,i])
  end
  return(haplotype)
end

function match(boolarray1::BitArray, boolarray2::BitArray, minimum_length::Int32)
  # TODO: Constrain the structure of the array - Static Types
  #Verify if there is match (size minimum_length) between two boolean sequences
  # TODO: Maybe you don't have to check the Whole array, depending on the size of the minimum match....
  # TODO: Also, maybe it's faster to detect if you start checking from the center of the sequence?
  maximum_match = 0
	size_of_sequence = size(boolarray2,1)
	current_maximum = 0
	for i = 1:size_of_sequence
		if boolarray1[i] != boolarray2[i]
			current_maximum+=1
			if maximum_match < current_maximum
				maximum_match = current_maximum
      end
		else
			current_maximum = 0
    end
  end

	return(maximum_match >= minimum_length)
end

function host_presents_peptide(peptide::BitArray, host::Host, minimum_length::Int32)
  # TODO: Constrain the structure of the array - Static Types
  for i = 1:2
    # TODO: Is there a better way to select the array?
    if match(host.sequences[:,i], peptide, minimum_length)
      return(true)
    end
  end
  return(false)
end

function to_Integer(boolean_array::BitArray)
  #Transforma um sequência binária num inteiro
  length = size(boolean_array, 1)
	integer = Int32(0)
	for i = 1:length
		if boolean_array[i]
			integer = integer + Int32(2^(i-1))
    end
  end
	return integer
end

function Heterosigozity(frequency::Array{Int32, 1})
  summation = sum(frequency)
  newfreq = frequency/summation
  return (1 - newfreq'*newfreq)[1]
end


# Main --------------------------------------------------------------------------

function main()

end

main()
