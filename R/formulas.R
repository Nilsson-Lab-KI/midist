#
# More documentation will follow
#


# Takes as input a broken formula and merges it back to a character string
merge_formula <- function(broken_formula) {
  broken_formula <- broken_formula[which(broken_formula$ElementNumber != 0), ]
  merged_formula <- as.vector(apply(as.data.frame(apply(broken_formula, 1, paste0, collapse = "")), 2, paste0, collapse = ""))
  return(gsub(" ", "", merged_formula))
}


#' Compares a formula to a list of formulas, and returns the match if there is any
compare_formula_diff_to_list <- function(formula_diff, allowed_reactions) {
  ind <- which(is.na(unlist(lapply(formula_diff, match, allowed_reactions))) != T)
  if (length(ind) > 0) {
    return(formula_diff[[ind]])
  }
}


#' Compares a mass to a list of masses, and returns the match if there is any within the given error range
compare_mass_diff_to_list <- function(mass_diff, allowed_masses, tol = 10) {
  error_l <- mass_diff - mass_diff * tol * 10^-6
  error_u <- mass_diff + mass_diff * tol * 10^-6

  ind <- which(allowed_masses > error_l & allowed_masses < error_u)
  if (length(ind) > 0) {
    return(allowed_masses[ind])
  }
}



# formulas is a list of length two, that stores the formulas to be compared to each other
# the output is a list that represents the reaction in both directions
get_formula_difference <- function(formulas) {
  broken_formulas <- lapply(formulas, break_formula)
  diff <- list()
  diff[[1]] <- data.frame(
    Symbol = broken_formulas[[1]]$Symbol,
    ElementNumber = broken_formulas[[1]]$ElementNumber - broken_formulas[[2]]$ElementNumber
  )
  diff[[2]] <- data.frame(
    Symbol = broken_formulas[[1]]$Symbol,
    ElementNumber = broken_formulas[[2]]$ElementNumber - broken_formulas[[1]]$ElementNumber
  )

  return(lapply(diff, merge_formula))
}

# get the periodic table elements (symbols)
Symbol <- c(
  "C", "H", "He", "Li", "Be", "B", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
  "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
  "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
  "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
  "Ts", "Og"
)
break_formula <- function(formula) {
  # make a data frame for the input formula
  brokenFormula <- as.data.frame(Symbol)
  brokenFormula$ElementNumber <- 0

  # splitting the formula into the smallest pieces
  formula_split <- unlist(strsplit(formula, ""))

  # merging the lower case letters (if there is any) to their upper case letter
  grep_lower <- grepl("[a-z]", formula_split)
  l <- which(grep_lower == TRUE)
  if (length(l) > 0) {
    formula_split[l - 1] <- paste(formula_split[l - 1], formula_split[l], sep = "")
    formula_split[l] <- NA
    formula_split <- formula_split[-which(is.na(formula_split) == TRUE)]
  }


  # find elements' location in the broken formula
  element_ind <- which(formula_split %in% Symbol)

  # now we need to assign "1" in case the last element doesn't have a number
  if (element_ind[length(element_ind)] == length(formula_split)) {
    formula_split <- c(formula_split, "1")
  }

  # now we need to assign "1" to any element in the middle that doesn't have a number
  if (length(element_ind) > 1) {
    for (i in 1:(length(element_ind) - 1)) {
      row_index <- which(brokenFormula$Symbol == formula_split[element_ind[i]])
      if (element_ind[i] + 1 == element_ind[i + 1]) {
        brokenFormula$ElementNumber[row_index] <- 1
      }
      if (element_ind[i] + 1 != element_ind[i + 1]) {
        brokenFormula$ElementNumber[row_index] <- as.numeric(paste(formula_split[(element_ind[i] + 1):(element_ind[i + 1] - 1)], collapse = ""))
      }
    }
  }

  i <- length(element_ind)
  row_index <- which(brokenFormula$Symbol == formula_split[element_ind[i]])
  brokenFormula$ElementNumber[row_index] <- as.numeric(paste(formula_split[(element_ind[i] + 1):length(formula_split)], collapse = ""))
  # brokenFormula <- brokenFormula[which(brokenFormula$ElementNumber!=0),]
  return(brokenFormula)
}


get_formula_mass <- function(formula, periodic_table) {
  broken_formula <- break_formula(formula)
  broken_formula$mass <- periodic_table$AtomicMass[match(broken_formula$Symbol, periodic_table$Symbol)] * broken_formula$ElementNumber
  return(sum(broken_formula$mass))
}
