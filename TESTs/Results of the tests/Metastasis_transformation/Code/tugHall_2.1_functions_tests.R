# Clone class
clone <- setRefClass(
    # name of the class
    Class = "Clone",
    # Field
    fields = list(
        id = "numeric",          # identificator
        parent = "numeric",      # parent ID (for first - 0)
        N_cells = "numeric",     # Number of cells in clone 
        c = "numeric",           # split counter
        d = "numeric",           # probability of division
        i = "numeric",           # probability of hayflick limit
        m = "numeric",           # probability that gene normal function is destroyed due to epigenome abnormality
        a = "numeric",           # probability of apoptosis
        s = "numeric",           # coefficient in the of apoptosis
        k = "numeric",           # probability of cell death by environment
        E = "numeric",           # coefficient of friction term against to the split probability,
                                 # coefficient for determination the max number of cells that can exist in the primary tumor (Nmax = 1/E)
        Nmax = "numeric",        # the max number of cells that can exist in the primary tumor (Nmax = 1/E)
        im = "numeric",          # invasion/ metastasis probability
        Ha = "numeric",          # apoptosis probability difference (apoptosis gene x weight) 
        Him = "numeric",         # invasion/ metastasis probability difference (invasion/ metastasis gene x weight) 
        Hi = "numeric",          # mitotic restriction probability difference (immortalized gene x weight)
        Hd = "numeric",          # divide probability difference (cancer gene x weight)
        Hb = "numeric",          # friction coefficient (angiogenic gene x weight)
        gene = "numeric",        # flag for cancer gene function deficit (for drivers)
        pasgene = "numeric",     # flag for cancer gene as passenger dysfunction 
        posdriver = "character", # position of cancer gene damage (function deficit)
        pospasngr = "character", # position of cancer gene damage (maintenance of function)
        mutden = "numeric",      # gene mutation density
        invasion = "logical",    # Wetting/Displacement flag:    TRUE: Wetting/Displacement      FALSE: Limited pattern
        birthday = "numeric"     # time step of birth of cell
    ),

    # Method
    methods = list(
        # Initialize
        initialize = function(gene_size, id=1, parent=0, c=0, d=d0, i=1, m=m0,  N_cells = 1, 
                              mutden=0, a=0, k=k0, E=E0, Nmax=0, gene=NULL, pasgene=NULL,
                              posdriver=NULL, pospasngr=NULL, invasion=FALSE, s=s0,birthday=0) {
            id <<- id
            parent <<- parent
            N_cells <<- N_cells
            c <<- c
            d <<- d
            i <<- i
            m <<- m
            s <<- s
            birthday <<- birthday
            mutden <<- mutden
            if (is.null(a)) {
                a <<- 1/(1+exp(-s*(mutden - 0.5)))
            } else {
                a <<- a
            }
            k <<- k
            E <<- E
            Nmax <<- 1.0 / E
            im <<- 0
            Ha <<- 0
            Him <<- 0
            Hi <<- 0
            Hd <<- 0
            Hb <<- 0

            if (is.null(gene)) {
                gene <<- rep(0, gene_size)
            } else {
                gene <<- gene
            }
            if (is.null(pasgene)) {
              pasgene <<- rep(0, gene_size)
            } else {
              pasgene <<- pasgene
            }
            if (is.null(posdriver)) {
                posdriver <<- rep("", gene_size)
            } else {
                posdriver <<- posdriver
            }
            if (is.null(pospasngr)) {
                pospasngr <<- rep("", gene_size)
            } else {
                pospasngr <<- pospasngr
            }
            invasion <<- invasion
        },
        # Apoptosis
        calcApoptosis = function() {
#            if (mutden <= sum(gene)/length(gene)) {
#                a1 = 1/(1+exp(s*(mutden - 0.5)))
#                mutden <<- sum(gene)/length(gene)
#                a2 = 1/(1+exp(s*(mutden - 0.5)))
#                a <<- a - (a1 - a2)

          mutden <<- sum(gene)/length(gene)
          a <<- 1/(1+exp(-1*s*(mutden - 0.5)))
                if (a < 0) {
                    a <<- 0
                }

        },
        # Aggregate
        calcMutden = function() {
            mutden <<- sum(gene)/length(gene)
        }
    )
)

# Environ class
environ <- setRefClass(
    # the class name
    Class = "Environ",

    # Fields
    fields = list(
        T = "numeric",           # time counter
        N = "numeric",           # localized clones number
        M = "numeric",           # number of infiltrting / metastatic clones
        F = "numeric",           # a coeffitient (Nmax = F/E) that determines the maximal number of cells 
                                 # that can exist in the primary tumor when the hallmark is engraved  
        c = "numeric",           # average number of divisions 
        d = "numeric",           # mean value of spliting probability
        i = "numeric",           # average value of spliting probability
        a = "numeric",           # average value of apoptosis probability
        k = "numeric",           # average probability of cell death
        E = "numeric",           # average value of coefficients of friction term proportional to N, for splitting probability
        Nmax = "numeric",        # Maximal number of cells that can exist 
        im = "numeric",          # average value of invasion / metastasis probability
        Ha = "numeric",
        Him = "numeric",
        Hi = "numeric",
        Hd = "numeric",
        Hb = "numeric",
        type = "numeric",        # invasion / metastatic ratio
        gene = "numeric",        # cancer gene damage rate
        posdriver = "character", # cancer gene damage position (function deficit)
        pospasngr = "character", # cancer gene damage position (maintaince of function)
        mutden = "numeric",      # average mutation rate
        last_id = "numeric"
    ),

    # Methods
    methods = list(
        # Initialize
        initialize = function(F0) {
            T <<- 0
            N <<- 0
            M <<- 0
            F <<- F0
        }
    )
)


# Cancer gene
oncogene <- setRefClass(
    # name of the class
    Class = "OncoGene",

    # Fields
    fields = list(
        name = "character",   # Cancer gene name list
        onsp = "character",   # oncogene/suppressor indicator
        cds = "numeric",      # cancer gene CDS base number list
        len = "numeric"       # number of cancer genes
    ),

    # Methods
    methods = list(
        # read the configuration file
        read = function(file) {
            data = read.table(file, sep="\t")
            name0 = NULL
            onsp0 = NULL
            cds0 = NULL
            for (i in 1:nrow(data)) {
                name <<- as.character(data[i, 1])
                if (!is.element(name, name0)) {
                    name0 = c(name0, name)
                    type = as.character(data[i, 4])
                    if (type == "?") {
                        if (runif(1) > 0.5) {
                            type = "o"
                        } else {
                            type = "s"
                        }
                    }
                    onsp0 = c(onsp0, type)
                    cds0 = c(cds0, as.numeric(as.character(data[i, 2])))
                }
            }
            name <<- name0
            onsp <<- onsp0
            cds <<- cds0
            len <<- length(name0)
        }
    )
)


hallmark <- setRefClass(
    # 
    Class = "HallMark",

    # 
    fields = list(
        Ha = "numeric",       # (Evading apoptosis)
        Hi = "numeric",       # (immortalization limit)
        Hd = "numeric",       # (Insensitivity to anti-growth signals || Self-sufficiency in growth signals)
        Hb = "numeric",       # (Sustained angiogenesis)
        Him = "numeric",      # (Tissue invasion & metastasis)
        Ha_w = "numeric",     # 
        Hi_w = "numeric",     # 
        Hd_w = "numeric",     # 
        Hb_w = "numeric",     # 
        Him_w = "numeric",    # 
        notHa = "numeric"
    ),

    # 
    methods = list(
        # 
        read = function(file, names) {
            data <- read.table(file, sep="\t")
            Ha0 = NULL
            Hi0 = NULL
            Hd0 = NULL
            Hb0 = NULL
            Him0 = NULL
            Ha0_w = NULL
            Hi0_w = NULL
            Hd0_w = NULL
            Hb0_w = NULL
            Him0_w = NULL
            w_flg = FALSE
            if (ncol(data) >= 5) {
                w_flg = TRUE
            }
            # Acquire gene name and Hallmark coefficient by function from definition file
            for (i in 1:nrow(data)) {
                if (data[i, 3] == "apoptosis") {
                    Ha0 = c(Ha0, as.character(data[i, 1]))
                    if (w_flg) {
                        Ha0_w = c(Ha0_w, as.numeric(as.character(data[i, 5])))
                    }
                } else if (data[i, 3] == "immortalization") {
                    Hi0 = c(Hi0, as.character(data[i, 1]))
                    if (w_flg) {
                        Hi0_w = c(Hi0_w, as.numeric(as.character(data[i, 5])))
                    }
                } else if (data[i, 3] == "anti-growth" | data[i, 3] == "growth") {
                    Hd0 = c(Hd0, as.character(data[i, 1]))
                    if (w_flg) {
                        Hd0_w = c(Hd0_w, as.numeric(as.character(data[i, 5])))
                    }
                } else if (data[i, 3] == "angiogenesis") {
                    Hb0 = c(Hb0, as.character(data[i, 1]))
                    if (w_flg) {
                        Hb0_w = c(Hb0_w, as.numeric(as.character(data[i, 5])))
                    }
                } else if (data[i, 3] == "invasion") {
                    Him0 = c(Him0, as.character(data[i, 1]))
                    if (w_flg) {
                        Him0_w = c(Him0_w, as.numeric(as.character(data[i, 5])))
                    }
                }
            }
            Ha <<- match(Ha0, names)
            notHa <<- setdiff(seq(1,length(names)),Ha)
            Hi <<- match(Hi0, names)
            Hd <<- match(Hd0, names)
            Hb <<- match(Hb0, names)
            Him <<- match(Him0, names)

            # if there is no Hallmark coefficient then generate a Hallmark coefficient as a random number - beta distribution
            if (!w_flg) {
                if (length(Ha) > 0) {
                    Ha_rnd = 1:length(Ha)
                } else {
                    Ha_rnd = c()
                }
                total0 = length(Ha)
                if (length(Hi) > 0) {
                    Hi_rnd = (total0 + 1):(total0 + length(Hi))
                } else {
                    Hi_rnd = c()
                }
                total0 = total0 + length(Hi)
                if (length(Hd) > 0) {
                    Hd_rnd = (total0 + 1):(total0 + length(Hd))
                } else {
                    Hd_rnd = c()
                }
                total0 = total0 + length(Hd)
                if (length(Hb) > 0) {
                    Hb_rnd = (total0 + 1):(total0 + length(Hb))
                } else {
                    Hb_rnd = c()
                }
                total0 = total0 + length(Hb)
                if (length(Him) > 0) {
                    Him_rnd = (total0 + 1):(total0 + length(Him))
                } else {
                    Him_rnd = c()
                }
                total = total0 + length(Him)
                # random = runif(total)
                random = rbeta(total, 0.01, 1)
                Ha0_w = random[Ha_rnd]
                Hi0_w = random[Hi_rnd]
                Hd0_w = random[Hd_rnd]
                Hb0_w = random[Hb_rnd]
                Him0_w = random[Him_rnd]
            }
            # Total by genetic mode 
            Ha_sum = sum(Ha0_w)
            Hi_sum = sum(Hi0_w)
            Hd_sum = sum(Hd0_w)
            Hb_sum = sum(Hb0_w)
            Him_sum = sum(Him0_w)
            Ha_w <<- Ha0_w/Ha_sum
            Hi_w <<- Hi0_w/Hi_sum
            Hd_w <<- Hd0_w/Hd_sum
            Hb_w <<- Hb0_w/Hb_sum
            Him_w <<- Him0_w/Him_sum
        },
        # Change the cell variables
        # mode = 2 Corresponding (Hallmark) Gene Mode
        updateClone = function(clone1, F) {
          # Apoptosis
          clone1$calcApoptosis()
          clone1$Ha = sum(clone1$gene[Ha]*Ha_w)
          clone1$a =  clone1$a - clone1$Ha
            if (clone1$a < 0) {
              clone1$a = 0
            }
          ### TEST to set a = 0 to exclude process 
          clone1$a = 0
          
            # Not dead - Immortalized
          clone1$Hi = sum(clone1$gene[Hi]*Hi_w)
          clone1$i = 1 - clone1$Hi
            if (clone1$i < 0) {
              clone1$i = 0
            }
            # Angiogenesis
          clone1$Hb = sum(clone1$gene[Hb]*Hb_w)
            
          clone1$E = E0 / (1 + F * clone1$Hb)
          clone1$Nmax = 1.0 / clone1$E
            
            # Cancer gene, tumor suppressor gene
          clone1$Hd = sum(clone1$gene[Hd]*Hd_w)
          clone1$Him = sum(clone1$gene[Him]*Him_w)
            
          clone1$d = clone1$Hd + d0     # d0 is defined in tugHall_2.1.R file 
                if (clone1$d > 1) {clone1$d = 1}
            if (!clone1$invasion) {
              clone1$d = clone1$d - clone1$E * env$N
              if (clone1$d < 0) {clone1$d = 0}
            }

                # Invasion metastasis
                if (!clone1$invasion) {
                  clone1$im = clone1$Him 
            } else {
              clone1$im = 1
            }
          ### Test for invasion* let im = 0.5
          clone1$im = 0.5
          
          },
        # Change the environment variables
        updateEnviron = function(env, clones) {
            sum_cell(env, clones)
            # env$M = sum()  # ceiling(length(clones) * env$type)
            # env$N = sum()  # length(clones)  - env$M
        }
    )
)

# Function to update Hallmark and variable after division or under initialization
update_Hallmarks <- function(clone1) {
    # Hallmark
    hall$updateClone(clone1, env$F)
}


# aggregate
sum_cell <- function(env, clones) {
    if (length(clones) > 0) {
        avg = apply(matrix(unlist(lapply(clones, sum_mutation)),ncol=length(clones)),1,sum)  #  /length(clones)
        env$c = avg[1] / (env$N + env$M)
        env$d = avg[2] / (env$N + env$M)
        env$i = avg[3] / (env$N + env$M)
        env$a = avg[4] / (env$N + env$M)
        env$k = avg[5] / (env$N + env$M)
        env$E = avg[6] / (env$N + env$M)
        env$Nmax = avg[7] / (env$N + env$M)
        env$im = avg[8] / (env$N + env$M)
        env$Ha = avg[9] / (env$N + env$M)
        env$Him = avg[10] / (env$N + env$M)
        env$Hi = avg[11] / (env$N + env$M)
        env$Hb = avg[12] / (env$N + env$M)
        env$Hd = avg[13] / (env$N + env$M)
        env$type = avg[14] / (env$N + env$M)
        env$mutden = avg[15] / (env$N + env$M)
    } else {
        env$M = 0
        env$N = 0
        env$c = 0
        env$d = 0
        env$i = 0
        env$a = 0
        env$k = 0
        env$E = 0
        env$Nmax = 0
        env$im = 0
        env$Ha = 0
        env$Him = 0
        env$Hi = 0
        env$Hb = 0
        env$Hd = 0
        env$type = 0
    }
    env$posdriver = rep("", length(onco$name))
    env$pospasngr = rep("", length(onco$name))
}


sum_mutation <- function(clone1) {
  return(c(clone1$c*clone1$N_cells,      clone1$d*clone1$N_cells,    clone1$i*clone1$N_cells, 
           clone1$a*clone1$N_cells,      clone1$k*clone1$N_cells,    clone1$E*clone1$N_cells,
           clone1$Nmax*clone1$N_cells,   clone1$im*clone1$N_cells,   clone1$Ha*clone1$N_cells, 
           clone1$Him*clone1$N_cells,    clone1$Hi*clone1$N_cells,   clone1$Hb*clone1$N_cells, 
           clone1$Hd*clone1$N_cells,     ifelse(clone1$invasion,1,0), 
           clone1$mutden*clone1$N_cells)) # , 
           # clone1$N_cells * ifelse(clone1$invasion,0,1),   clone1$N_cells * ifelse(clone1$invasion,1,0) ))
         
#           clone1$gene*clone1$N_cells) )
}

# To calculate N and M numbers - normal and metastasis cells
sum_N_M <- function(env, clones) {
    if (length(clones) > 0) {
        avg = apply(matrix(unlist(lapply(clones, number_N_M)),ncol=length(clones)),1,sum)  #  /length(clones)
        env$N = avg[1]
        env$M = avg[2]
        return(env$N + env$M)
    }
}

number_N_M <- function(clone1) {
    return( c(clone1$N_cells * ifelse(clone1$invasion,0,1),   clone1$N_cells * ifelse(clone1$invasion,1,0) ))
}



#####################################################################################################################################

# trial function
# The resul is a number of new cells, if N_New < 0 it means that the number is decreased.
trial <- function(clone1) { 
    
    # trial for Environmental death of cell 
    N_die = calc_binom(1, clone1$N_cells, clone1$k)   # The number of cells to die due to the Environmental death of cells in clone

    # Apoptosis trial
    N_die =  N_die + calc_binom(1, clone1$N_cells, clone1$a) 

    # invasion / metastasis trial 
    if (clone1$im > 0) {
        if (!clone1$invasion) {
            N_die = N_die + calc_binom(1, clone1$N_cells, (1 - clone1$im))
            clone1$invasion = ifelse(clone1$im == 1,TRUE,FALSE)
        }
    }

    # The new number of cells in the clone: 
    clone1$N_cells = ifelse( (clone1$N_cells - N_die) > 0, (clone1$N_cells - N_die) , 0)
    
    N_new = clone1$N_cells   # the initial number to split / before trial / - all cells "want" to split

    # Fragmentation restriction trial
    if (clone1$c > 50) {
        N_new = calc_binom(1, N_new, (1 - clone1$i))
    }
    
    # Devide trial
        
    N_new = calc_binom(1, N_new, clone1$d)   # The ! FINAL ! number of cells to split, this number is the output value of function "trial"

    if ( clone1$N_cells > 0 ) clone1$c = clone1$c + N_new / clone1$N_cells  # How to calculate the counter of division ? - As an average value of the counters for all cells !  

    clone1$N_cells = clone1$N_cells + N_new       # The number of cells are increased due to the splitting 
    
    N_new_clones = 0             # The number of new clones arising from clone1
    
    # p = clone1$m * sum(onco$cds) 
    # if (p < 1) N_new_clones = calc_binom(1, 2 * N_new, p ) else N_new_clones = 2 * N_new

    N_new_clones = calc_binom(1, 2 * N_new, 1 - p0 )   # coefficient 2, because of parent and children independent mutations 
    
    clone1$N_cells = ifelse( (clone1$N_cells - N_new_clones) > 0, (clone1$N_cells - N_new_clones) , 0)
    
    return(N_new_clones)
}



# mutagenesis trial
trial_mutagenesis <- function(clone1,num_mut) {
  
                            # num_mut is a number of mutations in this NEW clone
                            # length of onco - the number of genes 
    mut1 <- rep(0,onco$len)
  
    sm <- sample(index_gene,  num_mut,  replace = TRUE,  prob = prob_CDS)
    
    mut1[sm] <- 1
    
    if (sum(mut1) == 0) stop("The mutation is zero, that is incorrect", call. = TRUE)    # We have known that mutation must occur 
    
    posg = ceiling( runif(length(mut1[mut1==1])) * onco$cds[mut1==1] )    # The position of the site in the gene to mutate
    
    random2 = runif(length(mut1[mut1==1]))
    
    mut2 = ifelse((onco$onsp[mut1 == 1] == 'o' & uo >= random2) |
                    (onco$onsp[mut1 == 1] == 's' & us >= random2), 1, 0)
    
    clone1$gene[mut1 == 1] = ifelse(clone1$gene[mut1 == 1] == 1 | mut2 == 1, 1, 0)
    clone1$pasgene[mut1 == 1] = ifelse(clone1$pasgene[mut1 == 1] == 1 | mut2 == 0, 1, 0)
    
    # For drivers    
    clone1$posdriver[mut1 == 1][mut2 == 1] = ifelse(clone1$posdriver[mut1==1][mut2==1] == '',
                                                    paste(posg[mut2==1],env$T,sep = ":"), 
                                                    paste(clone1$posdriver[mut1==1][mut2==1], paste(posg[mut2==1],env$T,sep = ":"), sep=","))
    
    # For passangers
    clone1$pospasngr[mut1 == 1][mut2 == 0] = ifelse(clone1$pospasngr[mut1==1][mut2==0] == '',
                                                    paste(posg[mut2==0],env$T,sep = ":"), 
                                                    paste(clone1$pospasngr[mut1==1][mut2==0], paste(posg[mut2==0],env$T,sep = ":"), sep=","))
}



# to make one copy for clone1 in clone_init function
clone_copy <- function(clone1) {
    env$last_id = env$last_id + 1
    
    return(clone$new(id=env$last_id, parent=clone1$id, c=clone1$c, d=clone1$d, i=clone1$i, mutden=clone1$mutden,
                        a=clone1$a, k=clone1$k, E=clone1$E, Nmax=clone1$Nmax, gene=clone1$gene, pasgene=clone1$pasgene, 
                        posdriver=clone1$posdriver, pospasngr=clone1$pospasngr, 
                        invasion=clone1$invasion, s=clone1$s, birthday=env$T))
}



# write log file
write_log <- function(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, E0, F0, m, uo, us, s, k, 
                      censore_n, censore_t, d0) {
    data <- c("genefile", "clonefile", "geneoutfile", "cloneoutfile", "logoutfile",
              "E", "F", "m", "uo", "us", "s", "k", "censore_n", "censore_t", "d0")
    data <- rbind(data, c(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile,
                          E0, F0, m, uo, us, s, k, censore_n, censore_t, d0))
    write(data, logoutfile, ncolumn=2, sep="\t")
}


write_geneout <- function(outfile, hall) {
    data <- c(onco$name[hall$Ha], onco$name[hall$Hi], onco$name[hall$Hd], onco$name[hall$Hb], onco$name[hall$Him])
    data <- rbind(data, c(rep("apoptosis", length(onco$name[hall$Ha])),
                          rep("immortalization", length(onco$name[hall$Hi])),
                          rep("growth|anti-growth", length(onco$name[hall$Hd])),
                          rep("angiogenesis", length(onco$name[hall$Hb])),
                          rep("invasion", length(onco$name[hall$Him]))))
    data <- rbind(data, c(hall$Ha_w, hall$Hi_w, hall$Hd_w, hall$Hb_w, hall$Him_w))
    data <- rbind(data, c(onco$onsp[hall$Ha], onco$onsp[hall$Hi], onco$onsp[hall$Hd], onco$onsp[hall$Hb], 
                          onco$onsp[hall$Him]))
    write(data, outfile, ncolumn=4, sep="\t")
}


write_header_test <- function(outfile, env) {
  header <- c('Time', 'N_cells', 'Repeat_number', 'ID', 'ParentID:Birthday', 'c\'', 'd\'', 'i\'', 'im\'', 'a\'',
              'k\'', 'E\'', 'N', 'Nmax\'', 'M', 'Ha', 'Him', 'Hi', 'Hd', 'Hb', 'type', 'mut_den',
              paste("PosDriver:", onco$name, sep=""), paste("PosPasngr:", onco$name, sep="") )      #   , 'Clone number', 'Passengers Clone number', 'Mix Clone number')
  write(header, outfile, append=FALSE, ncolumn=length(header), sep="\t")
}

write_cloneout_test <- function(outfile, env, clones, isFirst, rp) {
  data <- c(env$T, '-', rp, '-', '-', env$c, env$d, env$i, env$im, env$a, env$k, env$E, env$N,
            env$Nmax, env$M, env$Ha, env$Him, env$Hi, env$Hd, env$Hb, env$type, env$mutden,
            env$gene,env$pospasngr, '-', '-', '-')
  write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
  
  # In the TESTs we need NOT to write each clone, only average values, that is why we use FALSE in the condition to skip this procedure
  if (FALSE & length(clones) > 0 & isFirst) {
    for (i in 1:length(clones)) {
      clone1 = clones[[i]]
      
      data <- c(env$T, clone1$N_cells, i, clone1$id, paste(clone1$parent,clone1$birthday,sep = ":"), clone1$c, clone1$d, 
                clone1$i, clone1$im, clone1$a, clone1$k, clone1$E, env$N, clone1$Nmax, env$M,
                clone1$Ha, clone1$Him, clone1$Hi, clone1$Hd, clone1$Hb, ifelse(clone1$invasion,1,0), 
                clone1$mutden, clone1$posdriver, clone1$pospasngr)           #   , vc, pasvc, mixvc)
      write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
    }
  }
}


write_weights <- function(outfile, hall) {
    #data <- c("Hallmarks", "Designation", onco$name)
    data <- data.frame( "Gene" = onco$name)   
    data$Gene <-   as.character(data$Gene)
    
    w <- rep(0.0, onco$len)
    for (j in 1:onco$len) { if ( length( which(j==hall$Ha) ) > 0  ) w[j] = hall$Ha_w[which(j==hall$Ha)]  }
    data <- cbind(data, "Apoptosis ($H_a$)" = w)
    
    w <- rep(0.0, onco$len)
    for (j in 1:onco$len) { if ( length( which(j==hall$Hb) ) > 0  ) w[j] = hall$Hb_w[which(j==hall$Hb)]  }
    data <- cbind(data, "Angiogenesis ($H_b$)" = w)
    
    w <- rep(0.0, onco$len)
    for (j in 1:onco$len) { if ( length( which(j==hall$Hd) ) > 0  ) w[j] = hall$Hd_w[which(j==hall$Hd)]  }
    data <- cbind(data, "Growth / Anti-growth ($H_d$)" = w)
    
    w <- rep(0.0, onco$len)
    for (j in 1:onco$len) { if ( length( which(j==hall$Hi) ) > 0  ) w[j] = hall$Hi_w[which(j==hall$Hi)]  }
    data <- cbind(data, "Immortalization ($H_i$)" = w)
    
    w <- rep(0.0, onco$len)
    for (j in 1:onco$len) { if ( length( which(j==hall$Him) ) > 0  ) w[j] = hall$Him_w[which(j==hall$Him)]  }
    data <- cbind(data, "Invasion / Metastasis ($H_{im}$)" = w)
    
    write.table(data, outfile, sep="\t", row.names = FALSE)
}


# initial clone setting
init_clones <- function(clonefile, clone1) {
    mpos <- regexpr("\\.", clonefile)[1]
    if (mpos != -1) {
        name <- substr(clonefile, 1, mpos - 1)
    } else {
        name <- clonefile
    }
    clones = NULL
    n <- as.numeric(name)
    if (!is.na(n) && is.numeric(n)) {
        factor = n / sum(clone1$m*onco$cds)
        f2 = 1.0
        while (TRUE) {
            if (sum(floor(clone1$m*onco$cds*factor*f2 + 0.5)) >= n) {
                break
            }
            f2 = f2 + 0.1
        }
        nums = floor(clone1$m*onco$cds*factor*f2 + 0.5)
        clones = NULL
        for (i in 1:n) {
            clones = c(clones, clone_copy(clone1))
        }
        pos = 0
        for (i in 1:length(nums)) {
            if (nums[i] > 0) {
                for (j in 1:nums[i]) {
                    if (pos + j <= n) {
                        clones[[pos + j]]$gene[i] = 1
                    }
                }
                pos = pos + nums[i]
            }
        }
    } else {
        data = read.table(clonefile, sep="\t")
        n <- nrow(data)
        
        for (i in 1:n) {
            clone2 = clone_copy(clone1)
            p <- match(onco$name, str_trim(strsplit(as.character(data[i,2]),",")[[1]]))
            clone2$gene[seq(1,length(onco$name))[!is.na(p)]] = 1
            clone2$N_cells = as.numeric(data[i,3])
            clones = c(clones, clone2)
        }
    }
    for (i in 1:n) {
        clones[[i]]$id = i
        clones[[i]]$parent = 0
        clones[[i]]$birthday = 0
        clones[[i]]$posdriver = ifelse(clones[[i]]$gene == 1,
                                      paste(ceiling(runif(onco$len)*onco$cds),"0",sep = ":"),
                                      clones[[i]]$posdriver)
        clones[[i]]$calcMutden()
        clones[[i]]$calcApoptosis()
    }
    env$last_id = n
    return(as.list(clones))
}

# Genetic recombination
read_w <- function(file) {
    if (!is.null(file) & !is.na(file)) {
        data <- read.table(file, sep="\t")
        w <- NULL
        r = 1
        start = 1
        while (TRUE) {
            window = ceiling(runif(1)*ncol(data))
            if (start + window > ncol(data)) {
                window = ncol(data) - start
            }
            w = c(w, data[r,start:(start+window)])
            if (start+window == ncol(data)) {
                break
            }
            start = start + window
            r = ifelse(r==1,2,1)                  
        }
        hall$setW(w)
    }
}


### Function to calculate binominal distribution including BIG NUMBERS like 10^12 and more using approximation with normal distribution 

calc_binom <- function(tr,n,p){
    if (n*p < 10^8){
        ou <- rbinom(tr,n,p) 
    } else {
        m <- n * p
        s <- sqrt(  n*p* (1-p)  )
        ou <- rnorm(tr,mean = m, sd = s)
    }
    
    return(  round( ou )  )
}


model_test <- function(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, 
                       E0, F0, m0, uo, us, s0, k0, censore_n, censore_t, d0, rp) {
  
    if (rp == 1) { write_log(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile, 
                           E0, F0, m0, uo, us, s0, k0, censore_n, censore_t, d0)  }  # write input parameters 
    onco = oncogene$new()        # make the vector onco about the hallmarks
    onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
    hall = hallmark$new()        # make a vector hall with hallmarks parameters
    hall$read(genefile, onco$name)     # read from the genefile - 'gene_cds2.txt'
    env = environ$new(F0)               # new vector for average values of cells
    assign("env", env, env=.GlobalEnv)
    assign("onco", onco, env=.GlobalEnv)
    assign("hall", hall, env=.GlobalEnv)
    assign("uo", uo, env=.GlobalEnv)
    assign("us", us, env=.GlobalEnv)
    clone1 = clone$new(gene_size=length(onco$cds),
                     m=m0, s=s0, k=k0, E=E0)          # clone1  -  empty object of clone
    clones = init_clones(clonefile, clone1)           # clones - the clones with hallmarks from cellfile - cellinit.txt - initial cells  
    if (rp == 1) { write_geneout(geneoutfile, hall)  }                # write the geneout.txt file with initial hallmarks 
    if (rp == 1) { write_weights("Output/Weights.txt", hall)    }             # write the weights of genes for hallmarks 
    if (rp == 1) { write_header_test(cloneoutfile, env)           }        # 
    
    cells_number <- sum_N_M(env, clones)                 # to calculate cells numbers - N,M 
    lapply(clones,update_Hallmarks)                     # to calculate the Hallmarks and probabilities for initial cells
    hall$updateEnviron(env, clones)                     # make averaging for cells 
    isFirst = TRUE
    write_cloneout_test(cloneoutfile, env, clones, isFirst, rp)      #  write initial clones

    # The initialization for mutation process
    
    prob_CDS <<- vector(mode = "numeric",length = onco$len)
    
    sumCDS <<- sum(onco$cds)
    
    prob_CDS <<- onco$cds / sumCDS
    
    p0 <<- (1 - m0) ^ sumCDS
    print( paste0("The probability of an absence of the mutations is p0 = ", as.character(p0) )) 
    
    #if (p0 >= 0 && p0 < 1) print( paste0("The probability of an absence of the mutations is p0 = ", as.character(p0) )  ) 
    #       else stop("The probability of an absence of the mutations is ", p0,  ", that is  incorrect!", call. = TRUE)
    
    index_gene <<- vector(mode = "numeric",length = onco$len)
    for (i in 1:onco$len) index_gene[i] <<- i
    
    while(length(clones) > 0 && censore_n > cells_number && env$T < censore_t) {

        k_old = length(clones)          # the number of clones from last step
        
        clones_new <- NULL 
        
        N_clones_new = unlist(lapply(clones, trial)) 
        
        survived_clones = NULL
        
        for (i in 1:k_old) {
            if (N_clones_new[i] > 0)  for (j in 1:N_clones_new[i])  clones_new = c(clones_new,clone_copy(clones[[i]]) ) 
            
            # To delete the clones with N_cells == 0, which are dead 
            if (clones[[i]]$N_cells == 0 )  survived_clones = c(survived_clones, FALSE)  else survived_clones = c(survived_clones, TRUE )   
        }
        
        
        # The number of mutations for each NEW clone
        N_new <- length(clones_new)
        
        if ( N_new > 0) {
            num_mut <- numeric(0)
            num_mut <- rztbinom(N_new, sumCDS, m0) # Numbers of mutations for new clones
            # To apply the mutagenesis only to new clones with a number of mutations: 
            for (nn in 1:N_new)  trial_mutagenesis( clones_new[[nn]], num_mut[nn]  )
        } 
        
        # the new generation = the survived clones + new_clones 
        clones = c(clones[survived_clones],clones_new)

        cells_number <- sum_N_M(env, clones)                 # to calculate cells numbers - N,M for next step
        lapply(clones,update_Hallmarks) 
        hall$updateEnviron(env, clones)                      # to average probabilities and hallmarks
        
        env$T = env$T + 1                                    # to next step
        
        write_cloneout_test(cloneoutfile, env, clones, isFirst, rp)
        #print(c(env$T,env$N,env$M,env$last_id, length(clones), "N_clones_new = ", N_clones_new))

    }
}



# Exchange the 3rd and 4th coloumns in the genefile 
changeCol <- function(genefile) {
    exchange = read.table(file =genefile, header = FALSE)
    Vec1 = exchange$V1
    Vec2 = exchange$V2
    Vec3 = exchange$V3
    Vec4 = exchange$V4
    Vec5 = exchange$V5
    exch_new = data.frame(Vec1,Vec2,Vec4,Vec3,Vec5)
    genefileNew = substr(genefile,1,4)
    write.table(exch_new, file = genefileNew, append = FALSE, quote = TRUE, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = FALSE, qmethod = c("escape", "double"),
                fileEncoding = "")
    
    return(genefileNew) 
}

