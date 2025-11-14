import math
import cmath

# ============================ funcoes basicas de matriz ============================

def eh_quase_zero(x, eps=1e-10):
    # verifica se o numero x eh muito proximo de zero (para tratar erros de arredondamento)
    return abs(x) < eps

def copiar_matriz(A):
    # cria uma copia independente da matriz a
    return [linha[:] for linha in A]

def matriz_identidade(n):
    # devolve a matriz identidade n x n
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

def dimensoes_matriz(A):
    # devolve o numero de linhas e colunas da matriz a
    return (len(A), len(A[0]) if A else 0)

def matriz_transposta(A):
    # devolve a transposta da matriz a
    m, n = dimensoes_matriz(A)
    return [[A[i][j] for i in range(m)] for j in range(n)]

def multiplicar_matrizes(A, B):
    # faz o produto c = a * b (se as dimensoes forem compativeis)
    m, n = dimensoes_matriz(A)
    n2, p = dimensoes_matriz(B)
    assert n == n2, "dimensoes incompativeis em multiplicar_matrizes"
    C = [[0.0 for _ in range(p)] for __ in range(m)]
    for i in range(m):
        for k in range(n):
            aik = A[i][k]
            if aik == 0:
                continue
            for j in range(p):
                C[i][j] += aik * B[k][j]
    return C

def multiplicar_matriz_vetor(A, v):
    # faz o produto entre matriz a e vetor coluna v
    m, n = dimensoes_matriz(A)
    assert len(v) == n, "dimensao do vetor nao compativel com a matriz"
    return [sum(A[i][j] * v[j] for j in range(n)) for i in range(m)]

def multiplicar_linha_por_escalar(M, i, escalar):
    # multiplica a linha i da matriz m por um escalar
    M[i] = [x * escalar for x in M[i]]

def trocar_linhas(M, i, k):
    # troca a linha i com a linha k na matriz m
    M[i], M[k] = M[k], M[i]

def somar_multiplo_de_linha(M, i, k, escalar):
    # faz linha_i = linha_i + escalar * linha_k
    M[i] = [M[i][j] + escalar * M[k][j] for j in range(len(M[i]))]

def aumentar_matriz(A, B):
    # monta a matriz aumentada [a | b] juntando as colunas
    m1, n1 = dimensoes_matriz(A)
    m2, n2 = dimensoes_matriz(B)
    assert m1 == m2, "matrizes devem ter o mesmo numero de linhas em aumentar_matriz"
    return [A[i] + B[i] for i in range(m1)]

def forma_escalonada_reduzida(M, eps=1e-10):
    # calcula a forma escalonada reduzida de m usando eliminacao de gauss jordan
    A = copiar_matriz(M)
    m, n = dimensoes_matriz(A)
    linha_pivo = 0
    colunas_pivo = []
    for j in range(n):
        # passo 1: procurar um bom pivo na coluna j
        indice_pivo = None
        maior_valor = 0.0
        for k in range(linha_pivo, m):
            if abs(A[k][j]) > maior_valor + eps:
                maior_valor = abs(A[k][j])
                indice_pivo = k
        if indice_pivo is None or eh_quase_zero(maior_valor, eps):
            # se nao achou pivo bom, pula a coluna
            continue
        # passo 2: trazer o pivo para a linha correta
        trocar_linhas(A, linha_pivo, indice_pivo)
        valor_pivo = A[linha_pivo][j]
        # passo 3: transformar o pivo em 1
        multiplicar_linha_por_escalar(A, linha_pivo, 1.0 / valor_pivo)
        # passo 4: zerar os outros elementos da coluna do pivo
        for k in range(m):
            if k != linha_pivo and not eh_quase_zero(A[k][j], eps):
                fator = -A[k][j]
                somar_multiplo_de_linha(A, k, linha_pivo, fator)
        colunas_pivo.append(j)
        linha_pivo += 1
        if linha_pivo == m:
            break
    return A, colunas_pivo

def posto_matriz(A, eps=1e-10):
    # conta quantas linhas nao nulas ha na forma escalonada reduzida (posto)
    R, _ = forma_escalonada_reduzida(A, eps)
    r = 0
    for linha in R:
        if any(abs(x) > eps for x in linha):
            r += 1
    return r

def inversa_matriz_quadrada(A, eps=1e-10):
    # calcula a inversa de uma matriz quadrada a usando gauss jordan
    n, n2 = dimensoes_matriz(A)
    assert n == n2, "matriz deve ser quadrada em inversa_matriz_quadrada"
    # passo 1: montar a matriz aumentada [a | i]
    AI = aumentar_matriz(copiar_matriz(A), matriz_identidade(n))
    i = 0
    # passo 2: aplicar gauss jordan no lado esquerdo ate virar identidade
    for j in range(n):
        indice_pivo = None
        for k in range(i, n):
            if abs(AI[k][j]) > eps:
                indice_pivo = k
                break
        if indice_pivo is None:
            raise ValueError("matriz singular, sem inversa")
        trocar_linhas(AI, i, indice_pivo)
        valor_pivo = AI[i][j]
        multiplicar_linha_por_escalar(AI, i, 1.0 / valor_pivo)
        for k in range(n):
            if k != i and abs(AI[k][j]) > eps:
                fator = -AI[k][j]
                somar_multiplo_de_linha(AI, k, i, fator)
        i += 1
    # passo 3: o lado direito agora eh a inversa
    inversa = [linha[n:] for linha in AI]
    return inversa

def matriz_menos_lambda_vezes_identidade(A, lam):
    # calcula a matriz a - lam * i
    n, m = dimensoes_matriz(A)
    assert n == m, "matriz_menos_lambda_vezes_identidade so funciona para matriz quadrada"
    B = copiar_matriz(A)
    for i in range(n):
        B[i][i] = B[i][i] - lam
    return B

def base_nucleo_matriz(A, eps=1e-10):
    # calcula uma base para o nucleo de a, ou seja, solucoes de a x = 0
    m, n = dimensoes_matriz(A)
    R, colunas_pivo = forma_escalonada_reduzida(A, eps)
    colunas_pivo = set(colunas_pivo)
    colunas_livres = [j for j in range(n) if j not in colunas_pivo]
    base = []
    # para cada variavel livre criamos um vetor da base do nucleo
    for coluna_livre in colunas_livres:
        v = [0.0] * n
        v[coluna_livre] = 1.0
        linha_atual = 0
        for coluna in range(n):
            # identifica qual coluna tem pivo em cada linha
            if linha_atual < m and eh_quase_zero(1 - R[linha_atual][coluna], eps) and coluna in colunas_pivo:
                # escreve a variavel basica em funcao das variaveis livres
                v[coluna] = -R[linha_atual][coluna_livre]
                linha_atual += 1
        base.append(v)
    if not base:
        # se nao ha variavel livre, o nucleo eh apenas o vetor zero
        base = [[0.0] * n]
    return base

def base_espaco_coluna(A, eps=1e-10):
    # calcula uma base para o espaco coluna com as colunas pivo da matriz a
    R, colunas_pivo = forma_escalonada_reduzida(A, eps)
    if not colunas_pivo:
        return []
    base = []
    for j in colunas_pivo:
        coluna = [A[i][j] for i in range(len(A))]
        base.append(coluna)
    return base

def imprimir_vetores(vetores, label=""):
    # imprime uma lista de vetores com um rotulo opcional
    if label:
        print(label)
    for v in vetores:
        print("   ", [float(x) if isinstance(x, (int, float)) else x for x in v])

# ============================ leitura de transformacoes ============================

def ler_transformacao_por_equacao():
    """le transformacao linear atraves das equacoes que definem T"""
    print("\n" + "="*60)
    print("  DEFINIR TRANSFORMACAO LINEAR POR EQUACAO")
    print("="*60)
    print("\nEscolha a dimensao da transformacao:")
    print("  1 - Transformacao 2D: T: R^2 -> R^2")
    print("  2 - Transformacao 3D: T: R^3 -> R^3")
    print("  3 - Dimensao customizada: T: R^n -> R^m")
    
    escolha = input("\nDigite sua escolha (1, 2 ou 3): ").strip()
    
    if escolha == "1":
        # 2D - T(x,y) = (a1*x + b1*y, a2*x + b2*y)
        print("\n" + "-"*60)
        print("  TRANSFORMACAO T: R^2 -> R^2")
        print("-"*60)
        print("\nVoce vai definir T(x,y) especificando cada componente da saida.")
        print("\nExemplo: se T(x,y) = (2x+3y, 4x-y), voce digita:")
        print("  Componente 1 de T(x,y): 2 3    <- coeficientes de x e y")
        print("  Componente 2 de T(x,y): 4 -1   <- coeficientes de x e y")
        print("\nAgora defina SUA transformacao T(x,y) = (?, ?):\n")
        A = []
        for i in range(2):
            while True:
                try:
                    linha_str = input(f"  Componente {i + 1} de T(x,y) - coeficientes [x, y]: ").strip().split()
                    linha = [float(x) for x in linha_str]
                    if len(linha) != 2:
                        print(f"  ERRO: esperava 2 coeficientes, recebeu {len(linha)}. Tente novamente.")
                        continue
                    A.append(linha)
                    break
                except ValueError:
                    print("  ERRO: digite apenas numeros (use ponto para decimal).")
        
        # Mostrar a equação resultante
        print(f"\nTransformacao definida:")
        print(f"  T(x,y) = ({A[0][0]}x + {A[0][1]}y, {A[1][0]}x + {A[1][1]}y)")
        return A
        
    elif escolha == "2":
        # 3D - T(x,y,z) = (a1*x + b1*y + c1*z, a2*x + b2*y + c2*z, a3*x + b3*y + c3*z)
        print("\n" + "-"*60)
        print("  TRANSFORMACAO T: R^3 -> R^3")
        print("-"*60)
        print("\nVoce vai definir T(x,y,z) especificando cada componente da saida.")
        print("\nExemplo: se T(x,y,z) = (x+2y, 3y-z, x+z), voce digita:")
        print("  Componente 1 de T(x,y,z): 1 2 0    <- coeficientes de x, y, z")
        print("  Componente 2 de T(x,y,z): 0 3 -1   <- coeficientes de x, y, z")
        print("  Componente 3 de T(x,y,z): 1 0 1    <- coeficientes de x, y, z")
        print("\nAgora defina SUA transformacao T(x,y,z) = (?, ?, ?):\n")
        A = []
        for i in range(3):
            while True:
                try:
                    linha_str = input(f"  Componente {i + 1} de T(x,y,z) - coeficientes [x, y, z]: ").strip().split()
                    linha = [float(x) for x in linha_str]
                    if len(linha) != 3:
                        print(f"  ERRO: esperava 3 coeficientes, recebeu {len(linha)}. Tente novamente.")
                        continue
                    A.append(linha)
                    break
                except ValueError:
                    print("  ERRO: digite apenas numeros (use ponto para decimal).")
        
        # Mostrar a equação resultante
        print(f"\nTransformacao definida:")
        print(f"  T(x,y,z) = ({A[0][0]}x + {A[0][1]}y + {A[0][2]}z,")
        print(f"              {A[1][0]}x + {A[1][1]}y + {A[1][2]}z,")
        print(f"              {A[2][0]}x + {A[2][1]}y + {A[2][2]}z)")
        return A
        
    elif escolha == "3":
        # Customizada - T: R^n -> R^m
        print("\n" + "-"*60)
        print("  TRANSFORMACAO T: R^n -> R^m (CUSTOMIZADA)")
        print("-"*60)
        n = int(input("\nDimensao do DOMINIO (n): ").strip())
        m = int(input("Dimensao do CONTRADOMINIO (m): ").strip())
        
        print(f"\nTransformacao T: R^{n} -> R^{m}")
        print(f"\nVoce vai definir T(x_1,...,x_{n}) especificando {m} componentes.")
        print(f"Cada componente eh uma combinacao linear de x_1, ..., x_{n}.\n")
        
        A = []
        for i in range(m):
            while True:
                try:
                    linha_str = input(f"  Componente {i + 1} - coeficientes de [x_1,...,x_{n}]: ").strip().split()
                    linha = [float(x) for x in linha_str]
                    if len(linha) != n:
                        print(f"  ERRO: esperava {n} coeficientes, recebeu {len(linha)}. Tente novamente.")
                        continue
                    A.append(linha)
                    break
                except ValueError:
                    print("  ERRO: digite apenas numeros (use ponto para decimal).")
        
        print(f"\nTransformacao T: R^{n} -> R^{m} definida com sucesso!")
        return A
        
    else:
        print("ERRO: Opcao invalida! Usando transformacao identidade 2x2 como padrao.")
        return [[1.0, 0.0], [0.0, 1.0]]

def ler_matriz_operador():
    """le matriz de operador linear T: R^n -> R^n (matriz quadrada)"""
    print("\n" + "="*60)
    print("  MATRIZ DO OPERADOR LINEAR T: R^n -> R^n")
    print("="*60)
    print("\nOBS: Operadores lineares precisam de MATRIZES QUADRADAS!")
    print("\nEscolha a dimensao do espaco:")
    print("  1 - Operador em R^2 (matriz 2x2)")
    print("  2 - Operador em R^3 (matriz 3x3)")
    
    escolha = input("\nDigite sua escolha (1 ou 2): ").strip()
    
    if escolha == "1":
        n = 2
    elif escolha == "2":
        n = 3
    else:
        print("ERRO: Opcao invalida! Usando matriz 2x2 como padrao.")
        n = 2
    
    print(f"\n" + "-"*60)
    print(f"  MATRIZ DO OPERADOR T: R^{n} -> R^{n}")
    print("-"*60)
    print(f"\nDigite a MATRIZ QUADRADA {n}x{n} que representa o operador T.")
    print("(Cada coluna representa como T age nos vetores base)")
    print("\nExemplo (matriz diagonal 2x2):")
    print("  linha 1: 2 0  <- T multiplica e_1 por 2")
    print("  linha 2: 0 3  <- T multiplica e_2 por 3")
    print("\nAgora digite SUA matriz do operador:\n")
    
    A = []
    for i in range(n):
        while True:
            try:
                linha_str = input(f"  Linha {i + 1} ({n} numeros separados por espaco): ").strip().split()
                linha = [float(x) for x in linha_str]
                if len(linha) != n:
                    print(f"  ERRO: esperava {n} numeros, recebeu {len(linha)}. Tente novamente.")
                    continue
                A.append(linha)
                break
            except ValueError:
                print("  ERRO: digite apenas numeros (use ponto para decimal).")
    
    print(f"\nMatriz {n}x{n} criada com sucesso!\n")
    return A

# ============================ tarefa 1 ============================

def produto_vetorial_3d(u, v):
    # calcula o produto vetorial em r3 entre u e v
    x1, y1, z1 = u
    x2, y2, z2 = v
    return [
        y1 * z2 - z1 * y2,
        z1 * x2 - x1 * z2,
        x1 * y2 - y1 * x2
    ]

def executar_tarefa_1():
    # tarefa 1:
    # w = {(x,y,z) em r3 | a x + b y + c z = 0}
    # objetivo: encontrar uma base e a dimensao de w
    print("\n" + "="*70)
    print("  TAREFA 1: SUBESPACO W = {(x,y,z) ∈ R³ | ax + by + cz = 0}")
    print("="*70)
    print("\nVoce vai definir um PLANO que passa pela origem em R³.")
    print("O plano sera descrito pela equacao: ax + by + cz = 0")
    print("\nExemplo: a=1, b=2, c=3 define o plano x + 2y + 3z = 0\n")
    
    # passo 1: ler os coeficientes a, b, c
    a = float(input("Digite o coeficiente 'a' (coeficiente de x): ").strip())
    b = float(input("Digite o coeficiente 'b' (coeficiente de y): ").strip())
    c = float(input("Digite o coeficiente 'c' (coeficiente de z): ").strip())
    vetor_normal = [a, b, c]

    print(f"\nEquacao do plano: {a}x + {b}y + {c}z = 0")
    print("-"*70)

    # passo 2: caso especial em que a equacao vira 0 = 0 e w = r3 inteiro
    if eh_quase_zero(a) and eh_quase_zero(b) and eh_quase_zero(c):
        print("\nCASO ESPECIAL: a = b = c = 0")
        print("   A equacao vira 0 = 0, logo W = R³ inteiro!")
        base = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        print("\nBase de W (base canonica de R³):")
        imprimir_vetores(base)
        print("\nDimensao de W: 3")
        return

    # passo 3: caso geral, w eh um plano passando pela origem (dimensao 2)

    # passo 3a: escolher um vetor u1 ortogonal ao vetor normal (a,b,c)
    if not (eh_quase_zero(a) and eh_quase_zero(b)):
        # se a ou b nao forem zero, u1 = (-b, a, 0) eh ortogonal a (a, b, c)
        u1 = [-b, a, 0.0]
    else:
        # se a = b = 0 e c != 0, o plano eh z = 0 e podemos usar u1 = (1,0,0)
        u1 = [1.0, 0.0, 0.0]

    # passo 3b: calcular u2 como produto vetorial entre o vetor normal e u1
    u2 = produto_vetorial_3d(vetor_normal, u1)

    # passo 3c: se u2 sair quase zero, trocamos u1 por outro vetor e recalculamos
    if eh_quase_zero(u2[0]) and eh_quase_zero(u2[1]) and eh_quase_zero(u2[2]):
        u1 = [0.0, 1.0, 0.0]
        u2 = produto_vetorial_3d(vetor_normal, u1)

    # passo 4: mostrar a base e a dimensao
    print("\nRESULTADOS:")
    print("   Base de W (dois vetores que geram o plano):")
    imprimir_vetores([u1, u2])
    print("\nDimensao de W: 2 (W eh um plano)")
    print("-"*70)

# ============================ tarefa 2 ============================

def executar_tarefa_2():
    # tarefa 2:
    # entrada: matriz a que representa uma transformacao linear t: r^n -> r^m
    # objetivos:
    #  - base do nucleo de t
    #  - dimensao do nucleo (nullidade)
    #  - base da imagem de t
    #  - dimensao da imagem (posto)
    print("\n" + "="*70)
    print("  TAREFA 2: NUCLEO E IMAGEM DE TRANSFORMACAO LINEAR")
    print("="*70)
    print("\nEsta tarefa calcula:")
    print("  • Nucleo de T (Ker(T)): vetores que T leva em zero")
    print("  • Imagem de T (Im(T)): vetores que sao resultado de T")
    print("  • Verifica o Teorema do Posto: dim(dominio) = posto + nulidade")
    
    # passo 1: ler a matriz da transformacao com escolha do tipo
    A = ler_transformacao_por_equacao()
    m, n = dimensoes_matriz(A)
    
    print(f"\nTransformacao T: R^{n} -> R^{m}")
    print("-"*70)

    # passo 2: calcular o nucleo resolvendo a x = 0
    nucleo = base_nucleo_matriz(A)

    # passo 3: calcular a imagem usando as colunas pivo (espaco coluna)
    espaco_coluna = base_espaco_coluna(A)
    posto = len(espaco_coluna)

    # passo 4: mostrar resultados
    print("\nRESULTADOS:")
    print("\n[1] NUCLEO DE T (solucoes de Ax = 0):")
    imprimir_vetores(nucleo, "   Base do Nucleo:")

    # trata o caso em que o nucleo eh apenas {0}
    if len(nucleo) == 1 and all(eh_quase_zero(x) for x in nucleo[0]):
        nullidade = 0
        print("   AVISO: Nucleo = {0} (apenas o vetor nulo)")
    else:
        nullidade = len(nucleo)

    print(f"\n   Dimensao do Nucleo (NULIDADE): {nullidade}")

    print("\n[2] IMAGEM DE T (espaco coluna da matriz):")
    imprimir_vetores(espaco_coluna, "   Base da Imagem:")
    print(f"\n   Dimensao da Imagem (POSTO): {posto}")

    # passo 5: conferir teorema do posto
    print("\n[3] VERIFICACAO DO TEOREMA DO POSTO:")
    print(f"   dim(dominio) = {n}")
    print(f"   posto + nulidade = {posto} + {nullidade} = {posto + nullidade}")
    if posto + nullidade == n:
        print("   Teorema verificado com sucesso!")
    else:
        print("   AVISO: pode haver erro numerico!")
    print("-"*70)

# ============================ tarefa 3 ============================

def ler_base_do_usuario(nome_base, tamanho_vetor, quantidade_vetores):
    # le do usuario uma base (conjunto de vetores) para um espaco
    print(f"\n" + "-"*60)
    print(f"  MATRIZ DE MUDANCA - BASE '{nome_base.upper()}'")
    print("-"*60)
    print(f"\nVoce vai digitar uma MATRIZ {quantidade_vetores}x{quantidade_vetores}")
    print(f"cujas COLUNAS sao os {quantidade_vetores} vetores da base {nome_base}.")
    print(f"\nCada vetor deve ter {tamanho_vetor} coordenadas.")
    print(f"\nExemplo para R2:")
    print(f"  Vetor 1: 1 0  <- primeira coluna = (1,0)")
    print(f"  Vetor 2: 0 1  <- segunda coluna = (0,1)")
    print(f"\nAgora digite OS VETORES da base {nome_base}:\n")
    
    lista_vetores = []
    for i in range(quantidade_vetores):
        linha_str = input(f"  Vetor {i + 1} da base {nome_base} ({tamanho_vetor} coordenadas): ").strip().split()
        v = [float(x) for x in linha_str]
        if len(v) != tamanho_vetor:
            raise ValueError(f"ERRO: esperava {tamanho_vetor} coordenadas, recebeu {len(v)}")
        lista_vetores.append(v)
    
    print(f"\nBase {nome_base} definida com sucesso!")
    # Construir matriz cujas colunas sao os vetores da base
    # lista_vetores[i] = vetor i (linha)
    # Precisamos de matriz onde coluna j = vetor j
    n = len(lista_vetores)
    m = len(lista_vetores[0])
    P = [[lista_vetores[j][i] for j in range(n)] for i in range(m)]
    return P

def executar_tarefa_3():
    # tarefa 3:
    # objetivo: achar a matriz de t em relacao a bases beta (dominio) e gama (contradominio)
    # formula: [t]_{gama <- beta} = p_gama^{-1} * a * p_beta
    print("\n" + "="*70)
    print("  TAREFA 3: MUDANCA DE BASE PARA TRANSFORMACAO LINEAR")
    print("="*70)
    print("\nEsta tarefa converte a matriz de T entre diferentes bases.")
    print("Formula: [T]_(γ←β) = P_γ^(-1) × A × P_β")
    print("\nOnde:")
    print("  • A: matriz de T na base canonica")
    print("  • β (beta): nova base do DOMINIO")
    print("  • γ (gama): nova base do CONTRADOMINIO")
    
    # passo 1: ler a matriz de t com escolha do tipo
    print("\n" + "="*70)
    print("[PASSO 1/3] MATRIZ A - TRANSFORMACAO NA BASE CANONICA")
    print("="*70)
    A = ler_transformacao_por_equacao()
    m, n = dimensoes_matriz(A)
    
    print(f"\nTransformacao T: R^{n} -> R^{m} (na base canonica)")
    print("-"*70)

    # passo 2: ler a base beta do dominio
    print("\n" + "="*70)
    print("[PASSO 2/3] BASE BETA (β) - NOVA BASE DO DOMINIO R^" + str(n))
    print("="*70)
    P_beta = ler_base_do_usuario("beta", n, n)      # matriz n x n

    # passo 3: ler a base gama do contradominio
    print("\n" + "="*70)
    print("[PASSO 3/3] BASE GAMA (γ) - NOVA BASE DO CONTRADOMINIO R^" + str(m))
    print("="*70)
    P_gama = ler_base_do_usuario("gama", m, m)      # matriz m x m

    # passo 4: verificar se gama realmente eh base (matriz invertivel)
    try:
        P_gama_inv = inversa_matriz_quadrada(P_gama)
    except Exception as e:
        print("\nERRO: Os vetores da base gama nao sao linearmente independentes!")
        print("   (A matriz nao e invertivel)")
        print(f"   Detalhes: {e}")
        return

    # passo 5: opcionalmente verificar se beta tambem eh base
    try:
        _ = inversa_matriz_quadrada(P_beta)
    except Exception as e:
        print("\nERRO: Os vetores da base beta nao sao linearmente independentes!")
        print("   (A matriz nao e invertivel)")
        print(f"   Detalhes: {e}")
        return

    # passo 6: aplicar a formula [t]_{gama <- beta} = p_gama^{-1} * a * p_beta
    temp = multiplicar_matrizes(A, P_beta)
    matriz_t_beta_gama = multiplicar_matrizes(P_gama_inv, temp)

    # passo 7: mostrar matriz final de t nas novas bases
    print("\nRESULTADO:")
    print("\nMatriz de T nas bases beta e gama:")
    print("   [T]_(gama<-beta) = P_gama^(-1) x A x P_beta")
    print()
    for linha in matriz_t_beta_gama:
        print("   ", [float(x) for x in linha])
    print("-"*70)

# ============================ tarefa 4 ============================

def autovalores_2x2(A):
    # calcula autovalores de uma matriz 2x2 usando traco e determinante
    a, b = A[0]
    c, d = A[1]
    traco = a + d
    determinante = a * d - b * c
    discriminante = traco * traco - 4 * determinante
    raiz_discriminante = cmath.sqrt(discriminante)
    l1 = (traco + raiz_discriminante) / 2
    l2 = (traco - raiz_discriminante) / 2
    return [l1, l2]

def decomposicao_qr_por_gram_schmidt(A):
    # faz decomposicao qr de uma matriz quadrada usando gram schmidt classico
    m, n = dimensoes_matriz(A)
    # aqui supomos que a eh m x m
    Q = [[0.0] * m for _ in range(m)]
    R = [[0.0] * m for _ in range(m)]

    # separa as colunas de a em lista v
    V = [[A[i][j] for i in range(m)] for j in range(m)]
    U = []  # vetores ortogonais

    for j in range(m):
        vj = V[j][:]  # copia da coluna j
        # passo 1: remover projetcoes em relacao aos vetores anteriores
        for k in range(j):
            uk = U[k]
            produto_interno = sum(vj[i] * uk[i] for i in range(m))
            R[k][j] = produto_interno
            for i in range(m):
                vj[i] -= produto_interno * uk[i]
        # passo 2: normalizar o vetor para ter norma 1
        norma_vj = math.sqrt(sum(vj[i] * vj[i] for i in range(m)))
        if eh_quase_zero(norma_vj):
            # se a norma for quase zero, devolve q = identidade e r = a (caso degenerado)
            return matriz_identidade(m), copiar_matriz(A)
        uj = [vj[i] / norma_vj for i in range(m)]
        U.append(uj)
        R[j][j] = norma_vj

    # passo 3: montar q usando os vetores u como colunas
    for i in range(m):
        for j in range(m):
            Q[i][j] = U[j][i]
    return Q, R

def autovalores_3x3_por_qr(A, iters=60):
    # aproxima autovalores de uma matriz 3x3 aplicando o algoritmo qr repetidas vezes
    Ak = copiar_matriz(A)
    for _ in range(iters):
        Q, R = decomposicao_qr_por_gram_schmidt(Ak)
        Ak = multiplicar_matrizes(R, Q)
    # no final, os autovalores aproximados aparecem na diagonal de ak
    return [Ak[0][0], Ak[1][1], Ak[2][2]]

def autovetores_para_autovalor(A, lam, eps=1e-8):
    # calcula autovetores associados a lam resolvendo (a - lam i) x = 0
    B = matriz_menos_lambda_vezes_identidade(A, lam)
    vetores = base_nucleo_matriz(B, eps)

    # remove vetor zero se aparecer
    vetores_limpos = []
    for v in vetores:
        if any(abs(x) > eps for x in v):
            vetores_limpos.append(v)
    if not vetores_limpos:
        # se nao achou nada, relaxa um pouco a tolerancia
        vetores = base_nucleo_matriz(B, eps * 10)
        vetores_limpos = [v for v in vetores if any(abs(x) > eps * 10 for x in v)]
    return vetores_limpos

def executar_tarefa_4():
    # tarefa 4:
    # entrada: matriz quadrada n x n (aqui usamos n = 2 ou 3)
    # objetivos:
    #  - mostrar a matriz do operador
    #  - calcular autovalores
    #  - calcular autoespacos e autovetores associados
    print("\n" + "="*70)
    print("  TAREFA 4: AUTOVALORES, AUTOVETORES E AUTOESPACOS")
    print("="*70)
    print("\nEsta tarefa calcula:")
    print("  • Autovalores λ: valores que satisfazem det(A - λI) = 0")
    print("  • Autovetores v: vetores que satisfazem Av = λv")
    print("  • Autoespacos: subespacos gerados pelos autovetores")
    
    # passo 1: ler a matriz do operador linear
    A = ler_matriz_operador()
    n, m = dimensoes_matriz(A)

    if n != m or n not in (2, 3):
        print("\nERRO: Este programa suporta apenas matrizes 2x2 ou 3x3.")
        return

    print("\nMatriz do operador T: R^{} -> R^{}:".format(n, n))
    print("-"*70)
    for linha in A:
        print("   ", [float(x) for x in linha])
    print("-"*70)

    # passo 2: calcular autovalores de acordo com o tamanho
    if n == 2:
        autovalores = autovalores_2x2(A)
    else:
        autovalores = autovalores_3x3_por_qr(A)

    print("\nAUTOVALORES ENCONTRADOS:")
    if n == 2:
        print("   (Calculados pela formula do traco e determinante)")
    else:
        print("   (Aproximados pelo algoritmo QR iterativo)")
    print()
    for i, lam in enumerate(autovalores, 1):
        print(f"   lambda_{i} = {lam}")

    # passo 3: para cada autovalor, calcular base do autoespaco
    print("\nAUTOESPACOS E AUTOVETORES:")
    print("   (Base do nucleo de A - lambda*I para cada autovalor)")
    print()
    autovalores_processados = []
    contador = 1
    for lam in autovalores:
        # agrupar autovalores muito proximos (especialmente no caso 3x3)
        ja_visto = False
        for lam_existente in autovalores_processados:
            if abs(lam - lam_existente) < 1e-6:
                ja_visto = True
                break
        if ja_visto:
            continue

        autovalores_processados.append(lam)
        vetores = autovetores_para_autovalor(A, lam)

        print(f"   [{contador}] Para lambda = {lam}:")
        if vetores:
            imprimir_vetores(vetores, "      Base do autoespaco:")
            print(f"      Dimensao do autoespaco: {len(vetores)}")
        else:
            print("      AVISO: Nao foi possivel encontrar autovetores")
            print("         (pode ser problema numerico ou autovalor complexo)")
        print()
        contador += 1
    
    print("-"*70)

# ============================ menu principal ============================

def menu_principal():
    # menu simples para escolher qual tarefa de algebra linear executar
    while True:
        print("\n" + "="*70)
        print("  PROGRAMA DE ALGEBRA LINEAR")
        print("="*70)
        print("\nEscolha uma tarefa:")
        print("  1 - Base e dimensao de subespaco W em R³")
        print("  2 - Nucleo e imagem de transformacao linear")
        print("  3 - Mudanca de base para transformacao linear")
        print("  4 - Autovalores, autovetores e autoespacos")
        print("  0 - Sair do programa")
        print("-"*70)
        
        opcao = input("\nDigite o numero da opcao desejada: ").strip()
        
        if opcao == "1":
            executar_tarefa_1()
        elif opcao == "2":
            executar_tarefa_2()
        elif opcao == "3":
            executar_tarefa_3()
        elif opcao == "4":
            executar_tarefa_4()
        elif opcao == "0":
            print("\n" + "="*70)
            print("  Programa encerrado. Ate mais!")
            print("="*70 + "\n")
            break
        else:
            print("\nERRO: Opcao invalida! Escolha um numero de 0 a 4.")

if __name__ == "__main__":
    menu_principal()
