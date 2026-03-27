
import re

class LoadFasta_FastaViewer:

    # Compilation for the ProteinId_dictionary and the lists with its values for each Protein
    def protein_dictionary(fastaFile, extra_pattern):
        """Gets fasta File Path and saves all informations about the proteins 
        in separate lists 

        Parameters
        ----------
        fastaFile : path to fasta file

        Returns
        -------
        for all proteins either Decoy or normal Lists/ dictionary.
        """
        dictKeyAccession = {}
        proteinList = []
        proteinNameList = []
        proteinOSList = []
        dictKeyAccessionDECOY = {}
        proteinListDECOY = []
        proteinNameListDECOY = []
        proteinOSListDECOY = []

        # if there is no accession, so we have no key
        default = 0

        #with open(fastaFile) as file_content:
        file_content = open(fastaFile)

        # Go through the fasta file, line by line
        nextline = file_content.readline()

        while(True):

            # if protein declaration (next protein) was found
            if nextline.startswith('>'):
                
                # safe information about the protein declaration for when inserting informations into lists
                protein_declaration = nextline

                check_if_reverse = protein_declaration.find('_rev')

                # check for extra pattern if used
                check_for_extra_pattern = -1
                if extra_pattern != '':
                    check_for_extra_pattern = protein_declaration.find(extra_pattern)

                # find upper and lower index of Protein Accession (ID)
                bounds = [m.start() for m in re.finditer(r'\|', protein_declaration)]

                # if a Protein Accesion (ID) was found
                if len(bounds) >= 1:

                    # if = 1, then the protein accession start at character 1 
                    if len(bounds) == 1:
                        key = (protein_declaration[1:bounds[0]])
                    else:
                        key = (protein_declaration[bounds[0]+1:bounds[1]])

                    descr_upper_index = protein_declaration.find('OS=')

                    os_upper_index = -1

                    # find out up to which index the OS goes
                    if descr_upper_index != -1:
                        # start searchin in seqs from index 'OS='
                        os_upper_index = protein_declaration.find('(', descr_upper_index)

                    os = 'not found'

                    # if no upper index for the OS was found use the line from "OS=" till end of line
                    if (os_upper_index == -1):
                        # if "OS=" was found, else stick with "os = not found"
                        if descr_upper_index != -1:
                            stringFrom_OS_TillEndOfLine = protein_declaration[descr_upper_index+3:]
                            listOfAllWordsTillEndOfLine = stringFrom_OS_TillEndOfLine.split()
                            os = listOfAllWordsTillEndOfLine[0] + ' ' + listOfAllWordsTillEndOfLine[1]

                    name = 'not found'

                    if descr_upper_index != -1:
                        if len(bounds) == 1:
                            name = (protein_declaration[bounds[0]+1:descr_upper_index])
                        else:
                            name = (protein_declaration[bounds[1]+1:descr_upper_index])
                        
                    protein_sequence_string = ""
                    nextline = file_content.readline()

                    # read file line by line, till new protein (begins with '>')
                    while not nextline.startswith('>'):
                        protein_sequence_string += nextline[:-1]
                        try:
                            nextline = next(file_content)
                        except Exception:
                            break

                    # set the values inside of the dictionary and the lists
                    # is decoy or reverse(!=-1) -> use seperate List to safe informations
                    if protein_declaration.startswith('>DECOY') or check_if_reverse != -1 or check_for_extra_pattern != -1:
                        dictKeyAccessionDECOY[key] = protein_sequence_string
                        proteinListDECOY.append(protein_sequence_string)
                        proteinNameListDECOY.append(name)

                        if (os_upper_index == -1 or descr_upper_index == -1):
                            proteinOSListDECOY.append(os)
                        else:
                            proteinOSListDECOY.append((protein_declaration[descr_upper_index+3:os_upper_index]))

                    # regular protein -> use regular list 
                    else:
                        dictKeyAccession[key] = protein_sequence_string
                        proteinList.append(protein_sequence_string)
                        proteinNameList.append(name)

                        if (os_upper_index == -1 or descr_upper_index == -1):
                            proteinOSList.append(os)
                        else:
                            proteinOSList.append((protein_declaration[descr_upper_index+3:os_upper_index]))
                    
                # if no Protein Accession (ID) was found
                else:
                    descr_upper_index = protein_declaration.find('OS=')

                    if descr_upper_index == -1:
                        key = str(default)
                        name = 'not found'
                        default = default + 1
                    else:    
                        key = protein_declaration.split("OS=")[0]
                        name = key
                    

                    os_upper_index = -1

                    # find out up to which index the OS goes
                    if descr_upper_index != -1:
                        # start searchin in seqs from index 'OS='
                        os_upper_index = protein_declaration.find('(', descr_upper_index)

                    os = 'not found'

                    # if no upper index for the OS was found use the line from "OS=" till end of line
                    if (os_upper_index == -1):
                        # if "OS=" was found, else stick with "os = not found"
                        if descr_upper_index != -1:
                            stringFrom_OS_TillEndOfLine = protein_declaration[descr_upper_index+3:]
                            listOfAllWordsTillEndOfLine = stringFrom_OS_TillEndOfLine.split()
                            os = listOfAllWordsTillEndOfLine[0] + ' ' + listOfAllWordsTillEndOfLine[1]

                    protein_sequence_string = ""
                    nextline = file_content.readline()

                    # read file line by line, till new protein (begins with '>')
                    while not nextline.startswith('>'):
                        protein_sequence_string += nextline[:-1]
                        try:
                            nextline = next(file_content)
                        except Exception:
                            break

                    # set the values inside of the dictionary and the lists
                    # is decoy -> use seperate List to safe informations
                    if protein_declaration.startswith('>DECOY') or check_if_reverse != -1 or check_for_extra_pattern != -1:
                        dictKeyAccessionDECOY[key] = protein_sequence_string
                        proteinListDECOY.append(protein_sequence_string)
                        proteinNameListDECOY.append(name)

                        if (os_upper_index == -1 or descr_upper_index == -1):
                            proteinOSListDECOY.append(os)
                        else:
                            proteinOSListDECOY.append((protein_declaration[descr_upper_index+3:os_upper_index]))

                    else:
                        dictKeyAccession[key] = protein_sequence_string
                        proteinList.append(protein_sequence_string)
                        proteinNameList.append(name)

                        if (os_upper_index == -1 or descr_upper_index == -1):
                            proteinOSList.append(os)
                        else:
                            proteinOSList.append((protein_declaration[descr_upper_index+3:os_upper_index]))
            
            else:
                break

        return dictKeyAccession, proteinList, proteinNameList, proteinOSList, dictKeyAccessionDECOY, proteinListDECOY, proteinNameListDECOY, proteinOSListDECOY         


        


def main():

    # for testing purposes only
    # use your own path for it
    dictKeyAccession, proteinList, proteinNameList, proteinOSList, dictKeyAccessionDECOY, proteinListDECOY, proteinNameListDECOY, proteinOSListDECOY = LoadFasta_FastaViewer.protein_dictionary(
        "path/to/file.fasta", "insert extra decoy pattern here")

    proteinnameSub = input("Name: ")

    for protein_name in proteinNameList:
        if proteinnameSub in protein_name:
            print("Proteinname: " + protein_name)

    # here the users enters the protein accession (key), or only a part of it
    # e.g.: 'P00761'
    protein_accession_maybe_sub_sequence = input("Bitte Protein accession angeben: ")

    # index starts with 0 (for the lists and the dictionary)
    # search for key
    for protein_accession in dictKeyAccession:
        if protein_accession_maybe_sub_sequence in protein_accession:
            index = list(dictKeyAccession).index(protein_accession)
            print("ID: " + list(dictKeyAccession.keys())[index])
            print("Protein: " + dictKeyAccession.get(protein_accession))
            print("Proteinname: " + proteinNameList[index])
            print("OS: " + proteinOSList[index])


if __name__ == "__main__":
    main()
