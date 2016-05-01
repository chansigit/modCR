# This file contains methods regarding file manipulation (data extraction)
import re

# Read and Quality pre-line markers (FASTQ)
READ_MARKER = "@"
QUAL_MARKER = "\+"

# Read entire file content to a list
def readFileContent(filename):
    # Create an handler to the file
    fileHandler = open(filename, 'r');    
    data = [];
    # read lines one by one
    line = fileHandler.readline();
    while(line != ''):
        # trim the line
        line = line.strip();
        # Add the line to the content
        data.append(line);
        # Read the next line
        line = fileHandler.readline();
    # Close handler
    fileHandler.close();
    return data;

# Read entire file content to a list
def readFileReadsThin(filename):
    # Create an handler to the file
    fileHandler = open(filename, 'r');
    # read lines one by one
    line = fileHandler.readline();
    while(line != ''):
        # trim the line
        line = line.strip();
        if not line.startswith(">"):
            # Add the line to the content
            yield line
        # Read the next line
        line = fileHandler.readline();
    # Close handler
    fileHandler.close();


# Read entire file content to a list
def readFileContentThin(filename):
    # Create an handler to the file
    fileHandler = open(filename, 'r');    
    # read lines one by one
    line = fileHandler.readline();
    while(line != ''):
        # trim the line
        line = line.strip();
        # Add the line to the content
        yield line
        # Read the next line
        line = fileHandler.readline();
    # Close handler
    fileHandler.close();

# Read a file of reads in a FASTQ format into data structures of reads and scores
def GetReadsData(readsFile, reads, scores):

    # Create an handler to the file
    fileHandler = open(readsFile, 'r');

    # Create regular expressions for read and quality markers
    readMarker      = re.compile(READ_MARKER);
    qualityMarker   = re.compile(QUAL_MARKER);
    
    # read lines and store reads and wualities
    line = fileHandler.readline();
    isRead = 0;
    isQuality = 0;
    
    while(line != ''):

        # trim the line
        line = line.strip();
        
        # Store line if read or quality score
        if(isRead):
            reads.append(line);
        elif(isQuality):
            scores.append(line);
            
        # Analyze if next line is a read
        if(readMarker.match(line) != None):
            isRead = 1;
        else:
            isRead = 0;

        # Analyze if next line is a quality score string
        if(qualityMarker.match(line) != None):
            isQuality = 1;
        else:
            isQuality = 0;
            
        # Read the next line
        line = fileHandler.readline();

    # Close handler
    fileHandler.close();

# Write a list into file (each element in a row)
def WriteListToFile(elements, fileName):

    # Create an handler to the file
    fileHandler = open(fileName, 'w');
    
    # Write every element in a row
    for element in elements:
        fileHandler.write(str(element) + "\n");
        
    # Close handler
    fileHandler.close();
   
# Write a list into file (each element in a row)
def AppendListToFile(elements, fileName):

    # Create an handler to the file
    fileHandler = open(fileName, 'a');
    
    # Write every element in a row
    for element in elements:
        fileHandler.write(element + "\n");
        
    # Close handler
    fileHandler.close();     
