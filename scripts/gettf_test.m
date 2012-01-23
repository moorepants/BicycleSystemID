u = 1:25;
y = 1:25;

na = 0
nb = 2
nk = 1

j = 0;
a = [];
b = [];
for i = (10 + (nk - 1)):(length(u) - (na + (nk - 1)) - 10),
    i
    j = j + 1
    a(j, 1:na) = y(i:i + na - 1)
    a(j, na + 1:na + nb) = u(i + na - nk - nb + 1:i + na - nk)
    b(j, 1) = y(i + na)
end

