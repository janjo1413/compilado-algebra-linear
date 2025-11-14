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
    print("\n=== tarefa 1: subespaco w = {(x,y,z) | a x + b y + c z = 0} ===")
    # passo 1: ler os coeficientes a, b, c
    a = float(input("digite o valor de a: ").strip())
    b = float(input("digite o valor de b: ").strip())
    c = float(input("digite o valor de c: ").strip())
    vetor_normal = [a, b, c]

    # passo 2: caso especial em que a equacao vira 0 = 0 e w = r3 inteiro
    if eh_quase_zero(a) and eh_quase_zero(b) and eh_quase_zero(c):
        print("caso especial: a = b = c = 0, entao w = r3 inteiro")
        base = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        print("uma base possivel para w eh a base canonica:")
        imprimir_vetores(base)
        print("dimensao de w: 3")
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
    print("uma base possivel para w eh formada pelos vetores:")
    imprimir_vetores([u1, u2])
    print("dimensao de w: 2")

# ============================ tarefa 2 ============================

def ler_matriz_do_usuario(mensagem="matriz (m x n)"):
    # le do usuario uma matriz de tamanho m x n
    print("\n" + mensagem)
    # passo 1: ler o tamanho da matriz
    m = int(input("numero de linhas m: ").strip())
    n = int(input("numero de colunas n: ").strip())
    A = []
    print("digite cada linha com os elementos separados por espaco (exemplo: 1 2 3)")
    # passo 2: ler linha por linha
    for i in range(m):
        linha_str = input(f"linha {i + 1}: ").strip().split()
        linha = [float(x) for x in linha_str]
        if len(linha) != n:
            raise ValueError("numero de colunas digitado nao confere com n")
        A.append(linha)
    return A

def executar_tarefa_2():
    # tarefa 2:
    # entrada: matriz a que representa uma transformacao linear t: r^n -> r^m
    # objetivos:
    #  - base do nucleo de t
    #  - dimensao do nucleo (nullidade)
    #  - base da imagem de t
    #  - dimensao da imagem (posto)
    print("\n=== tarefa 2: base, dimensao, nucleo e imagem de t ===")
    # passo 1: ler a matriz da transformacao
    A = ler_matriz_do_usuario("matriz de t: r^n -> r^m (cada coluna representa t(e_j))")
    m, n = dimensoes_matriz(A)

    # passo 2: calcular o nucleo resolvendo a x = 0
    nucleo = base_nucleo_matriz(A)

    # passo 3: calcular a imagem usando as colunas pivo (espaco coluna)
    espaco_coluna = base_espaco_coluna(A)
    posto = len(espaco_coluna)

    # passo 4: mostrar resultados
    print("\nbase do nucleo de t (solucoes de a x = 0):")
    imprimir_vetores(nucleo)

    # trata o caso em que o nucleo eh apenas {0}
    if len(nucleo) == 1 and all(eh_quase_zero(x) for x in nucleo[0]):
        nullidade = 0
    else:
        nullidade = len(nucleo)

    print(f"dimensao do nucleo de t (nullidade): {nullidade}")

    print("\nbase da imagem de t (espaco coluna da matriz):")
    imprimir_vetores(espaco_coluna)
    print(f"dimensao da imagem de t (posto da matriz): {posto}")

    # passo 5: conferir teorema do posto
    print(f"\nverificacao do teorema do posto: n = {n} e posto + nullidade = {posto} + {nullidade} = {posto + nullidade}")

# ============================ tarefa 3 ============================

def ler_base_do_usuario(nome_base, tamanho_vetor, quantidade_vetores):
    # le do usuario uma base (conjunto de vetores) para um espaco
    print(f"\nbase {nome_base}: informe {quantidade_vetores} vetores, cada um com {tamanho_vetor} coordenadas")
    lista_vetores = []
    for i in range(quantidade_vetores):
        linha_str = input(f"{nome_base} vetor {i + 1}: ").strip().split()
        v = [float(x) for x in linha_str]
        if len(v) != tamanho_vetor:
            raise ValueError("tamanho de vetor invalido na base")
        lista_vetores.append(v)
    # devolve matriz cujas colunas sao os vetores da base
    return matriz_transposta(lista_vetores)

def executar_tarefa_3():
    # tarefa 3:
    # objetivo: achar a matriz de t em relacao a bases beta (dominio) e gama (contradominio)
    # formula: [t]_{gama <- beta} = p_gama^{-1} * a * p_beta
    print("\n=== tarefa 3: matriz de t de beta para gama ===")
    print("entre com a matriz de t (m x n) nas bases canonicas")
    # passo 1: ler a matriz de t em base canonica
    A = ler_matriz_do_usuario("matriz a de t em base canonica")
    m, n = dimensoes_matriz(A)

    # passo 2: ler a base beta do dominio
    P_beta = ler_base_do_usuario("beta", n, n)      # matriz n x n

    # passo 3: ler a base gama do contradominio
    P_gama = ler_base_do_usuario("gama", m, m)      # matriz m x m

    # passo 4: verificar se gama realmente eh base (matriz invertivel)
    try:
        P_gama_inv = inversa_matriz_quadrada(P_gama)
    except Exception as e:
        print("erro: vetores informados para gama nao formam uma base (matriz nao invertivel):", e)
        return

    # passo 5: opcionalmente verificar se beta tambem eh base
    try:
        _ = inversa_matriz_quadrada(P_beta)
    except Exception as e:
        print("erro: vetores informados para beta nao formam uma base (matriz nao invertivel):", e)
        return

    # passo 6: aplicar a formula [t]_{gama <- beta} = p_gama^{-1} * a * p_beta
    temp = multiplicar_matrizes(A, P_beta)
    matriz_t_beta_gama = multiplicar_matrizes(P_gama_inv, temp)

    # passo 7: mostrar matriz final de t nas novas bases
    print("\nmatriz de t de beta para gama (p_gama^{-1} * a * p_beta):")
    for linha in matriz_t_beta_gama:
        print("   ", [float(x) for x in linha])

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
    print("\n=== tarefa 4: autovalores, autovetores e autoespacos ===")
    # passo 1: ler a matriz do operador linear
    A = ler_matriz_do_usuario("operador linear (matriz quadrada n x n, use n = 2 ou n = 3)")
    n, m = dimensoes_matriz(A)

    if n != m or n not in (2, 3):
        print("este programa suporta apenas matrizes 2x2 ou 3x3 na tarefa 4.")
        return

    print("\nmatriz do operador:")
    for linha in A:
        print("   ", [float(x) for x in linha])

    # passo 2: calcular autovalores de acordo com o tamanho
    if n == 2:
        autovalores = autovalores_2x2(A)
    else:
        autovalores = autovalores_3x3_por_qr(A)

    print("\nautovalores encontrados (em geral complexos para 2x2 e aproximados para 3x3):")
    for i, lam in enumerate(autovalores, 1):
        print(f"  lambda_{i} = {lam}")

    # passo 3: para cada autovalor, calcular base do autoespaco
    print("\nautoespacos e autovetores (base do nucleo de a - lambda i):")
    autovalores_processados = []
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

        print(f"\npara lambda = {lam}:")
        if vetores:
            imprimir_vetores(vetores, "   base do autoespaco:")
            print(f"   dimensao do autoespaco: {len(vetores)}")
        else:
            print("   nao foi possivel encontrar autovetores (pode ser problema numerico).")

# ============================ menu principal ============================

def menu_principal():
    # menu simples para escolher qual tarefa de algebra linear executar
    while True:
        print("\n==============================")
        print("menu algebra linear")
        print("==============================")
        print("1 - tarefa 1 (base e dimensao de w em r3)")
        print("2 - tarefa 2 (nucleo e imagem de t)")
        print("3 - tarefa 3 (matriz de t de beta para gama)")
        print("4 - tarefa 4 (autovalores e autovetores)")
        print("0 - sair")
        opcao = input("sua escolha: ").strip()
        if opcao == "1":
            executar_tarefa_1()
        elif opcao == "2":
            executar_tarefa_2()
        elif opcao == "3":
            executar_tarefa_3()
        elif opcao == "4":
            executar_tarefa_4()
        elif opcao == "0":
            print("encerrando o programa. ate mais.")
            break
        else:
            print("opcao invalida, tente novamente.")

if __name__ == "__main__":
    menu_principal()
