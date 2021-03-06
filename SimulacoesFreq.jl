include("ClassesFreq.jl")

# Parameters
Npep = Int32(20)
Nloci = Int32(1) #ALWAYS!
L = Int32(16) #tamanho das bitstrings
LT = Int32(7) #minimo de compatibilidade
Nhost = Int32(1000)
NG = Int32(10) #max number of pathogen genotypes
NS = Int32(5) #number os pathogen species
T = Int32(10000) #Geracoes de Hosts
muhost = 10.0^(-5)
mupath = 10.0^(-1)

# ------------------------------------------------------------------------------
#Misc
function initialize_random_host(length_of_locus::Int32)
  return(Host(bitrand(length_of_locus,2), 0))
end

function initialize_falses_host(length_of_locus::Int32)
  return(Host(falses(length_of_locus,2), 0))
end


function initialize_random_host(length_of_locus::Int32)
  return(Host(bitrand(length_of_locus,2), 0))
end

function initialize_random_host_population(number_of_hosts::Int32, mutation_rate::Float64, length_of_locus::Int32)
  host_population = fill(initialize_random_host(length_of_locus), number_of_hosts)
  for i = 2:number_of_hosts
    host_population[i] = initialize_random_host(length_of_locus)
  end
  return(PopHost(host_population, mutation_rate))
end

function initialize_falses_host_population(number_of_hosts::Int32, mutation_rate::Float64, length_of_locus::Int32)
  host_population = fill(initialize_falses_host(length_of_locus), number_of_hosts)
  for i = 2:number_of_hosts
    host_population[i] = initialize_falses_host(length_of_locus)
  end
  return(PopHost(host_population, mutation_rate))
end

function initialize_random_pathogen(length_of_locus::Int32, number_of_peptides::Int32)
  return(Path(bitrand(length_of_locus,number_of_peptides), 0))
end

function initialize_falses_path(length_of_locus::Int32, number_of_peptides::Int32)
  return(Path(falses(length_of_locus,number_of_peptides), 0))
end

function initialize_random_path_population(number_of_pathogens::Int32, mutation_rate::Float64, length_of_locus::Int32, number_of_peptides::Int32)
  random_number_of_pathogens = rand(1:number_of_pathogens)
  path_population = fill(initialize_random_pathogen(length_of_locus, number_of_peptides), random_number_of_pathogens)
  for i = 2:random_number_of_pathogens
    path_population[i] = initialize_random_pathogen(length_of_locus, number_of_peptides)
  end
  return(PopPath(path_population, mutation_rate))
end

function initialize_falses_path_population(number_of_pathogens::Int32, mutation_rate::Float64, length_of_locus::Int32, number_of_peptides::Int32)
  random_number_of_pathogens = rand(1:number_of_pathogens)
  path_population = fill(initialize_falses_path(length_of_locus, number_of_peptides), random_number_of_pathogens)
  for i = 2:random_number_of_pathogens
    path_population[i] = initialize_falses_path(length_of_locus, number_of_peptides)
  end
  return(PopPath(path_population, mutation_rate))
end

# ------------------------------------------------------------------------------
function main()
  for l = 1:50 #rodadas    #srand(17)
    println("Rodada = $(l)")

    #path_populations = fill(initialize_falses_path_population(NG, mupath, L, Npep), NS)
    #Calculating the mutation rates for each pathogen population
    #for i = 0:(NS-1)
  		#potencia = log10(mupath) + (log10(muhost) - log10(mupath))*i/(NS - 1)
  		#mu = 10.0^potencia
  		#path_populations[i+1] = initialize_falses_path_population(NG,mu,L,Npep)
    host_population = initialize_falses_host_population(Nhost, muhost, L)
    #  ----- Main loop -----
  	# For T generations
    #  ----- Main loop -----
  	# For T generations
    file = open("DRIFT_Rodada_=_$l.txt", "w")
  	H = fill(0.0,T)
    #counters = Array{Int32}[]
  	for i = 1:T
  		#println("T = $(i)")
  		 #For 10 generations of Paths
       #=
  		for j = 1:10

  			#INFECCAO
  			#For each host
  			for host in host_population.population_of_hosts
  				#contact with 1 pathoge from each spice
  				for k = 1:NS
  					pathogen = random_pathogen(path_populations[k], 1)
  					#verify if Path invades the Host (= compatibility in LT positions)
  					if pathogen_infects_host(pathogen[1],host,LT) #[1] porque a random_pathogen devolve um array de tamanho 1
  						#atualize o marcador do Path
  						infected!(pathogen[1])

  					else
  						#Host resisted
  						resisted!(host)
            end
          end
        end

  			#FITNESS, REPRODUCAO e MUTACAO - Paths
  			#For each PopPath
  			for k = 1:NS
  				path_populations[k] = reproduce(path_populations[k])
        end
      end
      =#

  		#FITNESS, REPRODUCAO e MUTACAO - Hosts
      H[i] = Heterosigozity(number_of_alleles(host_population))
  		#println(H[i])
      #println(fitness(host_population))
  		host_population = reproduce(host_population)
      write(file, "$(H[i]) ")
    end
    close(file)
  end
end
main()




#srand(17)
