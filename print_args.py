from sys import argv, stderr

print(f'\033[1;30margs \033[0m{len(argv)}\033[1;30m, argv [\033[0;33m', file=stderr, end='')
for i in range(len(argv)):
    print(f' {argv[i]}', file=stderr, end='')
print(" \033[1;30m]\033[0m", file=stderr)
