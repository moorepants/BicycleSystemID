function guesses = guess_nearby(seedGuess, numGuesses)
% function guesses = guess_nearby(seedGuess, numGuesses)
%
% Returns a series of random guesses normally distributed around a seed
% guess.
%
% Parameters
% ----------
% seedGuess : double, size(6)
%   A guess for the  six free parameters in the control system
%   identification.
% numGuesses : integer
%   The number of guesses to return.
%
% Returns
% -------
% guesses : double, size(numGuesses, 6)
%   A series of random guesses.

guesses = zeros(numGuesses, 6);
for i = 1:length(numGuesses)
    guesses(i, :) = normrnd(seedGuess, [50, 10, 10, 10, 10, 100]);
end
