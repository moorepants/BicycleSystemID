function guesses = guess_nearby(seedGuess, numGuesses)

guesses = zeros(numGuesses, 6);
for i = 1:length(numGuesses)
    guesses(i, :) = normrnd(seedGuess, [50, 10, 10, 10, 10, 100]);
end
