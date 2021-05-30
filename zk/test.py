import subprocess
import random

primes = [101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241]
for prime in primes:
    a = random.randint(1, prime-1)
    b = random.randint(1, prime-1)

    result = subprocess.run(['sage', 'prog.sage.py', str(a), str(b), str(prime)], stdout=subprocess.PIPE)
    myResult = result.stdout.decode('utf-8').strip()

    result = subprocess.run(['sage', 'truth.sage.py', str(a), str(b), str(prime)], stdout=subprocess.PIPE)
    otherResult = result.stdout.decode('utf-8').strip()

    print('[', a, ',', b, ',', prime, ']:', myResult, otherResult, myResult == otherResult)