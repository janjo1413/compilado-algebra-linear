import subprocess
import sys
import os

# caminho para o seu programa principal
# voce pode usar o caminho absoluto:
PROGRAM_PATH = r"C:\Users\jgque\OneDrive\Área de Trabalho\COMPILADOTRABALHOSALGEBRAJOAOGABRIELMIRANDAQUEIROZ\algebra_menu.py"
# se o script de teste estiver na mesma pasta do algebra_menu.py,
# voce pode simplificar para:
# PROGRAM_PATH = "algebra_menu.py"

PYTHON_EXE = sys.executable

# cada tupla: (nome_do_teste, entrada_que_sera_digitada_no_programa)
TESTS = [
    # ================= tarefa 1 =================
    # teste 1.1: a = b = c = 0 -> W = R^3
    (
        "tarefa1_w_R3",
        "1\n"     # escolhe opcao 1 no menu
        "0\n"     # a
        "0\n"     # b
        "0\n"     # c
        "0\n"     # volta ao menu e escolhe 0 para sair
    ),

    # teste 1.2: plano geral a=1, b=2, c=3
    (
        "tarefa1_plano_geral",
        "1\n"     # menu: tarefa 1
        "1\n"     # a
        "2\n"     # b
        "3\n"     # c
        "0\n"     # sair
    ),

    # ================= tarefa 2 =================
    # teste 2.1: identidade 2x2
    (
        "tarefa2_identidade_2x2",
        "2\n"         # menu: tarefa 2
        "2\n"         # m = 2
        "2\n"         # n = 2
        "1 0\n"       # linha 1
        "0 1\n"       # linha 2
        "0\n"         # sair
    ),

    # teste 2.2: matriz 2x2 com colunas dependentes
    (
        "tarefa2_dependente_2x2",
        "2\n"         # menu: tarefa 2
        "2\n"         # m = 2
        "2\n"         # n = 2
        "1 2\n"       # linha 1
        "2 4\n"       # linha 2
        "0\n"         # sair
    ),

    # teste 2.3: matriz 2x3, posto 2, nullidade 1
    (
        "tarefa2_2x3",
        "2\n"              # menu: tarefa 2
        "2\n"              # m = 2
        "3\n"              # n = 3
        "1 2 3\n"          # linha 1
        "4 5 6\n"          # linha 2
        "0\n"              # sair
    ),

    # ================= tarefa 3 =================
    # teste 3.1: bases canonicas (matriz nao muda)
    (
        "tarefa3_bases_canonicas",
        "3\n"              # menu: tarefa 3
        "2\n"              # m = 2
        "2\n"              # n = 2
        "1 2\n"            # linha 1 de A
        "3 4\n"            # linha 2 de A
        # base beta (2 vetores de 2 coords: base canonica)
        "1 0\n"            # beta vetor 1
        "0 1\n"            # beta vetor 2
        # base gama (2 vetores de 2 coords: base canonica)
        "1 0\n"            # gama vetor 1
        "0 1\n"            # gama vetor 2
        "0\n"              # sair
    ),

    # ================= tarefa 4 =================
    # teste 4.1: matriz 2x2 diagonal (autovalores 2 e 3)
    (
        "tarefa4_diagonal_2x2",
        "4\n"          # menu: tarefa 4
        "2\n"          # m = 2
        "2\n"          # n = 2
        "2 0\n"        # linha 1
        "0 3\n"        # linha 2
        "0\n"          # sair
    ),

    # teste 4.2: matriz 3x3 diagonal (autovalores 1,2,3)
    (
        "tarefa4_diagonal_3x3",
        "4\n"              # menu: tarefa 4
        "3\n"              # m = 3
        "3\n"              # n = 3
        "1 0 0\n"          # linha 1
        "0 2 0\n"          # linha 2
        "0 0 3\n"          # linha 3
        "0\n"              # sair
    ),
]


def rodar_um_teste(nome, entrada):
    """roda o programa algebra_menu.py com a entrada simulada."""
    print(f"rodando teste: {nome}...")
    resultado = subprocess.run(
        [PYTHON_EXE, PROGRAM_PATH],
        input=entrada,
        text=True,              # trabalha com strings (nao bytes)
        capture_output=True     # captura stdout e stderr
    )
    return resultado.stdout, resultado.stderr


def main():
    if not os.path.exists(PROGRAM_PATH):
        print("nao encontrei o arquivo do programa:")
        print("  ", PROGRAM_PATH)
        print("verifique se o caminho esta correto no script.")
        return

    with open("log_testes.txt", "w", encoding="utf-8") as log:
        log.write("LOG DE TESTES DO PROGRAMA algebra_menu.py\n")
        log.write("=========================================\n\n")

        for nome, entrada in TESTS:
            log.write("=" * 60 + "\n")
            log.write(f"TESTE: {nome}\n")
            log.write("-" * 60 + "\n")
            log.write("entrada enviada ao programa (simulando digitacao):\n")
            log.write(repr(entrada) + "\n\n")

            stdout, stderr = rodar_um_teste(nome, entrada)

            log.write("saida completa do programa (stdout):\n")
            log.write(stdout)
            log.write("\n")

            if stderr:
                log.write("mensagens de erro (stderr):\n")
                log.write(stderr)
                log.write("\n")

            log.write("\n\n")

    print("\npronto! arquivo 'log_testes.txt' foi gerado na mesma pasta do script.")
    print("voce pode abrir esse arquivo e me mandar o conteudo aqui para eu conferir.")


if __name__ == "__main__":
    main()
