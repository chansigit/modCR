# Longest Increasing subsequence Algorithm
def LIS(sequence):
   # Variable indicating the length of the longest subsequence
   optimal = 0;
   # List that maintain for every index J the position K of the smallest value 
   # in the sequence such that there is an increasing subsequence of length J 
   # ending at SEQ[K] on the range K (Note: J<=K<=Current Index)
   optimalEnds = [-1] * (len(sequence) + 1);
   # List maintaing the predecessor of the element in the longest subsequence 
   # ending at that position
   preds = [-1] * len(sequence);

   # Go over sequence indices
   for index in xrange(len(sequence)):
      # Search for the largest sequence ending before current element
      longestSubSoFar = LongestSunsequenceSoFar(sequence, optimalEnds, index);
      # Update the predecessor of current index
      preds[index] = optimalEnds[longestSubSoFar];
      # Update optimal ends and length of optimal subsequnce (if needed)
      if ((longestSubSoFar == optimal) or (sequence[index] < sequence[optimalEnds[longestSubSoFar + 1]])):
         optimalEnds[longestSubSoFar + 1] = index;
         optimal = max(optimal, longestSubSoFar + 1);
        
   sub = [];
   indices = [];
   if (optimal > 0):
      current = optimalEnds[optimal];
      sub.append(sequence[current]);
      indices.append(current);
      while (preds[current] != -1):
         current = preds[current];
         sub.append(sequence[current]);
         indices.append(current);
   sub.reverse(); 
   indices.reverse();
   return indices;
     
# Get the length of the longest subsequence that can end up to current index
def LongestSunsequenceSoFar(sequence, optimalEnds, index):
   # Current longest subsequnce
   result = 0;
   # Get the current number
   currentNumber = sequence[index];
   # Define the initial search boundaries
   low = 1;
   high = index;
   # Go over the array
   while (low <= high):
      current = (low + high) / 2;
      # If there is a valid subsequence of length current - search longer one
      endPosition = optimalEnds[current];
      if ((endPosition > -1) and (sequence[endPosition] < currentNumber)):
         low = current + 1;
         result = current;
      else:
         high = current - 1;
   return result;
 
   
#l = [1,2,3,4,5];
#l = [1,2,10,4,5];
#l = [100, 101, 1, 2, 3];
#LIS(l);
     
     

  