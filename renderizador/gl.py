#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# pylint: disable=invalid-name

"""
Biblioteca Gráfica / Graphics Library.

Desenvolvido por: <SEU NOME AQUI>
Disciplina: Computação Gráfica
Data: <DATA DE INÍCIO DA IMPLEMENTAÇÃO>
"""

import time         # Para operações com tempo
import gpu          # Simula os recursos de uma GPU
import math         # Funções matemáticas
import numpy as np  # Biblioteca do Numpy

class GL:
    """Classe que representa a biblioteca gráfica (Graphics Library)."""

    width = 800   # largura da tela
    height = 600  # altura da tela
    near = 0.01   # plano de corte próximo
    far = 1000    # plano de corte distante

    @staticmethod
    def setup(width, height, near=0.01, far=1000):
        """Definr parametros para câmera de razão de aspecto, plano próximo e distante."""
        GL.width = width
        GL.height = height
        GL.near = near
        GL.far = far

    @staticmethod
    def polypoint2D(point, colors):
        """Função usada para renderizar Polypoint2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Polypoint2D
        # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
        # de pontos x, y sempre na ordem. Assim point[0] é o valor da coordenada x do
        # primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é a
        # coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
        # pelo tamanho da lista e assuma que sempre vira uma quantidade par de valores.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polypoint2D
        # você pode assumir inicialmente o desenho dos pontos com a cor emissiva (emissiveColor).

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Polypoint2D : pontos = {0}".format(point)) # imprime no terminal pontos
        print("Polypoint2D : colors = {0}".format(colors)) # imprime no terminal as cores

        # Implementação simples: desenhar cada ponto
        emissive_color = colors.get('emissiveColor', [1.0, 1.0, 1.0])
        r = int(emissive_color[0] * 255)
        g = int(emissive_color[1] * 255)
        b = int(emissive_color[2] * 255)
        
        for i in range(0, len(point), 2):
            x = int(point[i])
            y = int(point[i + 1])
            if 0 <= x < GL.width and 0 <= y < GL.height:
                gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, [r, g, b])
        
    @staticmethod
    def polyline2D(lineSegments, colors):
        """Função usada para renderizar Polyline2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Polyline2D
        # Nessa função você receberá os pontos de uma linha no parâmetro lineSegments, esses
        # pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o valor da
        # coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto. Já point[2] é
        # a coordenada x do segundo ponto e assim por diante. Assuma a quantidade de pontos
        # pelo tamanho da lista. A quantidade mínima de pontos são 2 (4 valores), porém a
        # função pode receber mais pontos para desenhar vários segmentos. Assuma que sempre
        # vira uma quantidade par de valores.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Polyline2D
        # você pode assumir inicialmente o desenho das linhas com a cor emissiva (emissiveColor).

        print("Polyline2D : lineSegments = {0}".format(lineSegments)) # imprime no terminal
        print("Polyline2D : colors = {0}".format(colors)) # imprime no terminal as cores
        
        # Implementação simples: desenhar linhas conectando pontos consecutivos
        emissive_color = colors.get('emissiveColor', [1.0, 1.0, 1.0])
        r = int(emissive_color[0] * 255)
        g = int(emissive_color[1] * 255)
        b = int(emissive_color[2] * 255)
        
        for i in range(0, len(lineSegments) - 2, 2):
            x0 = int(lineSegments[i])
            y0 = int(lineSegments[i + 1])
            x1 = int(lineSegments[i + 2])
            y1 = int(lineSegments[i + 3])
            
            # Desenhar linha simples entre dois pontos
            GL._draw_line_simple(x0, y0, x1, y1, [r, g, b])

    @staticmethod
    def circle2D(radius, colors):
        """Função usada para renderizar Circle2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#Circle2D
        # Nessa função você receberá um valor de raio e deverá desenhar o contorno de
        # um círculo.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o Circle2D
        # você pode assumir o desenho das linhas com a cor emissiva (emissiveColor).

        print("Circle2D : radius = {0}".format(radius)) # imprime no terminal
        print("Circle2D : colors = {0}".format(colors)) # imprime no terminal as cores
        
        # Exemplo:
        pos_x = GL.width//2
        pos_y = GL.height//2
        gpu.GPU.draw_pixel([pos_x, pos_y], gpu.GPU.RGB8, [255, 0, 255])  # altera pixel (u, v, tipo, r, g, b)
        # cuidado com as cores, o X3D especifica de (0,1) e o Framebuffer de (0,255)


    @staticmethod
    def triangleSet2D(vertices, colors):
        """Função usada para renderizar TriangleSet2D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry2D.html#TriangleSet2D
        # Nessa função você receberá os vertices de um triângulo no parâmetro vertices,
        # esses pontos são uma lista de pontos x, y sempre na ordem. Assim point[0] é o
        # valor da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto.
        # Já point[2] é a coordenada x do segundo ponto e assim por diante. Assuma que a
        # quantidade de pontos é sempre multiplo de 3, ou seja, 6 valores ou 12 valores, etc.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, para o TriangleSet2D
        # você pode assumir inicialmente o desenho das linhas com a cor emissiva (emissiveColor).
        print("TriangleSet2D : vertices = {0}".format(vertices)) # imprime no terminal
        print("TriangleSet2D : colors = {0}".format(colors)) # imprime no terminal as cores

        # Implementação simples: desenhar triângulos preenchidos
        emissive_color = colors.get('emissiveColor', [1.0, 1.0, 1.0])
        r = int(emissive_color[0] * 255)
        g = int(emissive_color[1] * 255)
        b = int(emissive_color[2] * 255)
        
        for i in range(0, len(vertices), 6):
            if i + 5 < len(vertices):
                x1 = int(vertices[i])
                y1 = int(vertices[i + 1])
                x2 = int(vertices[i + 2])
                y2 = int(vertices[i + 3])
                x3 = int(vertices[i + 4])
                y3 = int(vertices[i + 5])
                
                # Preencher triângulo simples
                GL._fill_triangle_simple(x1, y1, x2, y2, x3, y3, [r, g, b])


    @staticmethod
    def triangleSet(point, colors):
        """Função usada para renderizar TriangleSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#TriangleSet
        # Nessa função você receberá pontos no parâmetro point, esses pontos são uma lista
        # de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x do
        # primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e
        # assim por diante.
        # No TriangleSet os triângulos são informados individualmente, assim os três
        # primeiros pontos definem um triângulo, os três próximos pontos definem um novo
        # triângulo, e assim por diante.
        # O parâmetro colors é um dicionário com os tipos cores possíveis, você pode assumir
        # inicialmente, para o TriangleSet, o desenho das linhas com a cor emissiva
        # (emissiveColor), conforme implementar novos materias você deverá suportar outros
        # tipos de cores.

        if not point or len(point) < 9:
            return

        # cor emissiva (X3D em 0..1 → framebuffer 0..255)
        rgb = colors.get('emissiveColor', [1.0, 1.0, 1.0]) if isinstance(colors, dict) else [1.0, 1.0, 1.0]
        col = [int(max(0, min(255, round(c*255)))) for c in rgb]

        # garante matrizes padrão caso ainda não tenham sido definidas
        M = getattr(GL, "_M", np.eye(4, dtype=float))
        V = getattr(GL, "_V", np.eye(4, dtype=float))
        P = getattr(GL, "_P", np.eye(4, dtype=float))
        PVM = P @ V @ M

        w, h = GL.width, GL.height

        def ndc_to_screen(ndc_xy):
            x_ndc, y_ndc = ndc_xy[0], ndc_xy[1]
            x = (x_ndc + 1.0) * 0.5 * (w - 1)
            y = (1 - (y_ndc + 1.0) * 0.5) * (h - 1)  # origem (0,0) no topo
            return int(round(x)), int(round(y))

        for i in range(0, len(point), 9):
            if i + 8 >= len(point):  # segurança
                break

            v0 = np.array([point[i],   point[i+1], point[i+2], 1.0], dtype=float)
            v1 = np.array([point[i+3], point[i+4], point[i+5], 1.0], dtype=float)
            v2 = np.array([point[i+6], point[i+7], point[i+8], 1.0], dtype=float)

            c0 = PVM @ v0
            c1 = PVM @ v1
            c2 = PVM @ v2

            # divisão por w -> NDC
            if c0[3] == 0 or c1[3] == 0 or c2[3] == 0:
                continue
            ndc0 = c0[:3] / c0[3]
            ndc1 = c1[:3] / c1[3]
            ndc2 = c2[:3] / c2[3]

            # mapeia para tela
            x0, y0 = ndc_to_screen(ndc0)
            x1, y1 = ndc_to_screen(ndc1)
            x2, y2 = ndc_to_screen(ndc2)

            # preenche o triângulo usando seu raster simples
            GL._fill_triangle_simple(x0, y0, x1, y1, x2, y2, col)

        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        # gpu.GPU.draw_pixel([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

    @staticmethod
    def viewpoint(position, orientation, fieldOfView):
        """Função usada para renderizar (na verdade coletar os dados) de Viewpoint."""
        # Na função de viewpoint você receberá a posição, orientação e campo de visão da
        # câmera virtual. Use esses dados para poder calcular e criar a matriz de projeção
        # perspectiva para poder aplicar nos pontos dos objetos geométricos.

        def T(tx, ty, tz):
            M = np.eye(4, dtype=float)
            M[:3, 3] = [tx, ty, tz]
            return M

        def R_axis_angle(x, y, z, angle):
            v = np.array([x, y, z], dtype=float)
            n = np.linalg.norm(v)
            if n == 0.0:
                return np.eye(4, dtype=float)
            x, y, z = v / n
            c, s, C = math.cos(angle), math.sin(angle), 1.0 - math.cos(angle)
            R3 = np.array([
                [x*x*C + c,   x*y*C - z*s, x*z*C + y*s],
                [y*x*C + z*s, y*y*C + c,   y*z*C - x*s],
                [z*x*C - y*s, z*y*C + x*s, z*z*C + c  ]
            ], dtype=float)
            R = np.eye(4, dtype=float)
            R[:3, :3] = R3
            return R
        # ----------------------

        px, py, pz = (position if position else [0.0, 0.0, 10.0])
        if orientation and len(orientation) == 4:
            ox, oy, oz, ang = orientation
        else:
            ox, oy, oz, ang = 0.0, 1.0, 0.0, 0.0

        # inverso da rotação = rotacionar por -ang no mesmo eixo
        R_inv = R_axis_angle(ox, oy, oz, -ang)
        T_inv = T(-px, -py, -pz)

        # View matrix (leva mundo -> câmera)
        V = R_inv @ T_inv
        setattr(GL, "_V", V)

        # Projeção perspectiva a partir do FOV vertical
        width, height = GL.width, GL.height
        aspect = width / max(1, height)
        fovy = float(fieldOfView) if fieldOfView else math.radians(60.0)
        f = 1.0 / math.tan(fovy * 0.5)

        z_near = GL.near
        z_far  = GL.far

        P = np.array([
            [f/aspect, 0,  0,                                0],
            [0,        f,  0,                                0],
            [0,        0,  (z_far+z_near)/(z_near-z_far),   (2*z_far*z_near)/(z_near-z_far)],
            [0,        0, -1,                                0]
        ], dtype=float)

        setattr(GL, "_P", P)
    @staticmethod
    def transform_in(translation, scale, rotation):
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_in será chamada quando se entrar em um nó X3D do tipo Transform
        # do grafo de cena. Os valores passados são a escala em um vetor [x, y, z]
        # indicando a escala em cada direção, a translação [x, y, z] nas respectivas
        # coordenadas e finalmente a rotação por [x, y, z, t] sendo definida pela rotação
        # do objeto ao redor do eixo x, y, z por t radianos, seguindo a regra da mão direita.
        # ESSES NÃO SÃO OS VALORES DE QUATÉRNIOS AS CONTAS AINDA PRECISAM SER FEITAS.
        # Quando se entrar em um nó transform se deverá salvar a matriz de transformação dos
        # modelos do mundo para depois potencialmente usar em outras chamadas. 
        # Quando começar a usar Transforms dentre de outros Transforms, mais a frente no curso
        # Você precisará usar alguma estrutura de dados pilha para organizar as matrizes.

        def T(tx, ty, tz):
            M = np.eye(4, dtype=float)
            M[:3, 3] = [tx, ty, tz]
            return M

        def S(sx, sy, sz):
            return np.diag([sx, sy, sz, 1.0]).astype(float)

        def R_axis_angle(x, y, z, angle):
            v = np.array([x, y, z], dtype=float)
            n = np.linalg.norm(v)
            if n == 0.0:
                return np.eye(4, dtype=float)
            x, y, z = v / n
            c, s, C = math.cos(angle), math.sin(angle), 1.0 - math.cos(angle)
            R3 = np.array([
                [x*x*C + c,   x*y*C - z*s, x*z*C + y*s],
                [y*x*C + z*s, y*y*C + c,   y*z*C - x*s],
                [z*x*C - y*s, z*y*C + x*s, z*z*C + c  ]
            ], dtype=float)
            R = np.eye(4, dtype=float)
            R[:3, :3] = R3
            return R
        # -----------------------------

        tx, ty, tz = (translation if translation else [0.0, 0.0, 0.0])
        sx, sy, sz = (scale       if scale       else [1.0, 1.0, 1.0])
        if rotation and len(rotation) == 4:
            rx, ry, rz, ang = rotation
        else:
            rx, ry, rz, ang = 0.0, 0.0, 1.0, 0.0

        local = T(tx, ty, tz) @ R_axis_angle(rx, ry, rz, ang) @ S(sx, sy, sz)

        # inicializa pilha e M atual se ainda não existirem
        if not hasattr(GL, "_M"):
            setattr(GL, "_M", np.eye(4, dtype=float))
        if not hasattr(GL, "_M_stack"):
            setattr(GL, "_M_stack", [])

        # empilha M corrente e compõe com a transformação local
        GL._M_stack.append(GL._M.copy())
        GL._M = GL._M @ local  # compõe para Transform dentro de Transform

    @staticmethod
    def transform_out():
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_out será chamada quando se sair em um nó X3D do tipo Transform do
        # grafo de cena. Não são passados valores, porém quando se sai de um nó transform se
        # deverá recuperar a matriz de transformação dos modelos do mundo da estrutura de
        # pilha implementada.

        # restaura a matriz do topo da pilha
        if hasattr(GL, "_M_stack") and len(GL._M_stack) > 0:
            GL._M = GL._M_stack.pop()
        else:
            # se algo sair do esperado, volta para identidade
            GL._M = np.eye(4, dtype=float)

    @staticmethod
    def triangleStripSet(point, stripCount, colors):
        """Função usada para renderizar TriangleStripSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#TriangleStripSet
        # A função triangleStripSet é usada para desenhar tiras de triângulos interconectados,
        # você receberá as coordenadas dos pontos no parâmetro point, esses pontos são uma
        # lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor da coordenada x
        # do primeiro ponto, point[1] o valor y do primeiro ponto, point[2] o valor z da
        # coordenada z do primeiro ponto. Já point[3] é a coordenada x do segundo ponto e assim
        # por diante. No TriangleStripSet a quantidade de vértices a serem usados é informado
        # em uma lista chamada stripCount (perceba que é uma lista). Ligue os vértices na ordem,
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        # Cor emissiva como cor sólida do strip
        rgb = colors.get('emissiveColor', [1.0, 1.0, 1.0]) if isinstance(colors, dict) else [1.0, 1.0, 1.0]
        solid = [int(max(0, min(255, round(c*255)))) for c in rgb]

        # Helpers locais
        w, h = GL.width, GL.height
        M = getattr(GL, "_M", np.eye(4))
        V = getattr(GL, "_V", np.eye(4))
        P = getattr(GL, "_P", np.eye(4))
        PVM = P @ V @ M

        def to_screen(v3):
            v = np.array([v3[0], v3[1], v3[2], 1.0], dtype=float)
            clip = PVM @ v
            if clip[3] == 0: 
                return None
            ndc = clip[:3] / clip[3]
            x = (ndc[0] + 1.0) * 0.5 * (w - 1)
            y = (1.0 - (ndc[1] + 1.0) * 0.5) * (h - 1)  # origem no topo
            return int(round(x)), int(round(y))

        def draw_tri(i0, i1, i2):
            v0 = point[3*i0:3*i0+3]
            v1 = point[3*i1:3*i1+3]
            v2 = point[3*i2:3*i2+3]
            p0 = to_screen(v0); p1 = to_screen(v1); p2 = to_screen(v2)
            if p0 and p1 and p2:
                GL._fill_triangle_simple(p0[0], p0[1], p1[0], p1[1], p2[0], p2[1], solid)

        # Varre as tiras
        base = 0
        for count in stripCount:
            # vértices da tira: base .. base+count-1
            for k in range(count - 2):
                i0 = base + k
                i1 = base + k + 1
                i2 = base + k + 2
                # alterna orientação para manter winding consistente
                if k % 2 == 0:
                    draw_tri(i0, i1, i2)
                else:
                    draw_tri(i1, i0, i2)
            base += count

    @staticmethod
    def indexedTriangleStripSet(point, index, colors):
        """Função usada para renderizar IndexedTriangleStripSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/rendering.html#IndexedTriangleStripSet
        # A função indexedTriangleStripSet é usada para desenhar tiras de triângulos
        # interconectados, você receberá as coordenadas dos pontos no parâmetro point, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim point[0] é o valor
        # da coordenada x do primeiro ponto, point[1] o valor y do primeiro ponto, point[2]
        # o valor z da coordenada z do primeiro ponto. Já point[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedTriangleStripSet uma lista informando
        # como conectar os vértices é informada em index, o valor -1 indica que a lista
        # acabou. A ordem de conexão será de 3 em 3 pulando um índice. Por exemplo: o
        # primeiro triângulo será com os vértices 0, 1 e 2, depois serão os vértices 1, 2 e 3,
        # depois 2, 3 e 4, e assim por diante. Cuidado com a orientação dos vértices, ou seja,
        # todos no sentido horário ou todos no sentido anti-horário, conforme especificado.

        # Cor sólida
        rgb = colors.get('emissiveColor', [1.0, 1.0, 1.0]) if isinstance(colors, dict) else [1.0, 1.0, 1.0]
        solid = [int(max(0, min(255, round(c*255)))) for c in rgb]

        # Helpers
        w, h = GL.width, GL.height
        M = getattr(GL, "_M", np.eye(4))
        V = getattr(GL, "_V", np.eye(4))
        P = getattr(GL, "_P", np.eye(4))
        PVM = P @ V @ M

        def to_screen_by_idx(ii):
            v3 = point[3*ii:3*ii+3]
            v = np.array([v3[0], v3[1], v3[2], 1.0], dtype=float)
            clip = PVM @ v
            if clip[3] == 0:
                return None
            ndc = clip[:3] / clip[3]
            x = (ndc[0] + 1.0) * 0.5 * (w - 1)
            y = (1.0 - (ndc[1] + 1.0) * 0.5) * (h - 1)
            return int(round(x)), int(round(y))

        def draw_tri(i0, i1, i2):
            p0 = to_screen_by_idx(i0); p1 = to_screen_by_idx(i1); p2 = to_screen_by_idx(i2)
            if p0 and p1 and p2:
                GL._fill_triangle_simple(p0[0], p0[1], p1[0], p1[1], p2[0], p2[1], solid)

        # Processa blocos separados por -1
        strip = []
        for idx in index + [-1]:
            if idx == -1:
                if len(strip) >= 3:
                    for k in range(len(strip) - 2):
                        i0, i1, i2 = strip[k], strip[k+1], strip[k+2]
                        if k % 2 == 0:
                            draw_tri(i0, i1, i2)
                        else:
                            draw_tri(i1, i0, i2)
                strip = []
            else:
                strip.append(idx)

    @staticmethod
    def indexedFaceSet(coord, coordIndex, colorPerVertex, color, colorIndex,
                       texCoord, texCoordIndex, colors, current_texture):
        """Função usada para renderizar IndexedFaceSet."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#IndexedFaceSet
        # A função indexedFaceSet é usada para desenhar malhas de triângulos. Ela funciona de
        # forma muito simular a IndexedTriangleStripSet porém com mais recursos.
        # Você receberá as coordenadas dos pontos no parâmetro cord, esses
        # pontos são uma lista de pontos x, y, e z sempre na ordem. Assim coord[0] é o valor
        # da coordenada x do primeiro ponto, coord[1] o valor y do primeiro ponto, coord[2]
        # o valor z da coordenada z do primeiro ponto. Já coord[3] é a coordenada x do
        # segundo ponto e assim por diante. No IndexedFaceSet uma lista de vértices é informada
        # em coordIndex, o valor -1 indica que a lista acabou.
        # A ordem de conexão não possui uma ordem oficial, mas em geral se o primeiro ponto com os dois
        # seguintes e depois este mesmo primeiro ponto com o terçeiro e quarto ponto. Por exemplo: numa
        # sequencia 0, 1, 2, 3, 4, -1 o primeiro triângulo será com os vértices 0, 1 e 2, depois serão
        # os vértices 0, 2 e 3, e depois 0, 3 e 4, e assim por diante, até chegar no final da lista.
        # Adicionalmente essa implementação do IndexedFace aceita cores por vértices, assim
        # se a flag colorPerVertex estiver habilitada, os vértices também possuirão cores
        # que servem para definir a cor interna dos poligonos, para isso faça um cálculo
        # baricêntrico de que cor deverá ter aquela posição. Da mesma forma se pode definir uma
        # textura para o poligono, para isso, use as coordenadas de textura e depois aplique a
        # cor da textura conforme a posição do mapeamento. Dentro da classe GPU já está
        # implementadado um método para a leitura de imagens.

        # Os prints abaixo são só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("IndexedFaceSet : ")
        if coord:
            print("\tpontos(x, y, z) = {0}, coordIndex = {1}".format(coord, coordIndex))
        print("colorPerVertex = {0}".format(colorPerVertex))
        if colorPerVertex and color and colorIndex:
            print("\tcores(r, g, b) = {0}, colorIndex = {1}".format(color, colorIndex))
        if texCoord and texCoordIndex:
            print("\tpontos(u, v) = {0}, texCoordIndex = {1}".format(texCoord, texCoordIndex))
        if current_texture:
            image = gpu.GPU.load_texture(current_texture[0])
            print("\t Matriz com image = {0}".format(image))
            print("\t Dimensões da image = {0}".format(image.shape))
        print("IndexedFaceSet : colors = {0}".format(colors))  # imprime no terminal as cores

        # ---------- Helpers de pipeline ----------
        w, h = GL.width, GL.height
        M = getattr(GL, "_M", np.eye(4))
        V = getattr(GL, "_V", np.eye(4))
        P = getattr(GL, "_P", np.eye(4))
        PVM = P @ V @ M

        def to_screen_idx(ii):
            v3 = coord[3*ii:3*ii+3]
            v = np.array([v3[0], v3[1], v3[2], 1.0], dtype=float)
            clip = PVM @ v
            if clip[3] == 0:
                return None, None
            ndc = clip[:3] / clip[3]
            x = (ndc[0] + 1.0) * 0.5 * (w - 1)
            y = (1.0 - (ndc[1] + 1.0) * 0.5) * (h - 1)
            return (int(round(x)), int(round(y))), ndc[2]  # z_ndc, se precisar no futuro

        # Raster sólido (uma cor)
        em = colors.get('emissiveColor', [1.0, 1.0, 1.0]) if isinstance(colors, dict) else [1.0, 1.0, 1.0]
        solid = [int(max(0, min(255, round(c*255)))) for c in em]

        # Raster com cor por vértice (baricêntricas simples)
        def fill_triangle_color(p0, p1, p2, c0, c1, c2):
            minx = max(0, min(p0[0], p1[0], p2[0]))
            maxx = min(w-1, max(p0[0], p1[0], p2[0]))
            miny = max(0, min(p0[1], p1[1], p2[1]))
            maxy = min(h-1, max(p0[1], p1[1], p2[1]))

            def edge(a, b, c):
                return (c[0]-a[0])*(b[1]-a[1]) - (c[1]-a[1])*(b[0]-a[0])

            area = edge(p0, p1, p2)
            if area == 0:
                return
            for ypix in range(miny, maxy+1):
                for xpix in range(minx, maxx+1):
                    P = (xpix + 0.5, ypix + 0.5)
                    w0 = edge(p1, p2, P)
                    w1 = edge(p2, p0, P)
                    w2 = edge(p0, p1, P)
                    if (area > 0 and w0 >= 0 and w1 >= 0 and w2 >= 0) or \
                       (area < 0 and w0 <= 0 and w1 <= 0 and w2 <= 0):
                        a0 = w0/area; a1 = w1/area; a2 = w2/area
                        r = int(max(0, min(255, round(a0*c0[0] + a1*c1[0] + a2*c2[0]))))
                        g = int(max(0, min(255, round(a0*c0[1] + a1*c1[1] + a2*c2[1]))))
                        b = int(max(0, min(255, round(a0*c0[2] + a1*c1[2] + a2*c2[2]))))
                        gpu.GPU.draw_pixel([xpix, ypix], gpu.GPU.RGB8, [r, g, b])

        # ---------- Triangulação em fan ----------
        face = []
        for idx in coordIndex + [-1]:
            if idx == -1:
                if len(face) >= 3:
                    v0 = face[0]
                    for j in range(1, len(face)-1):
                        v1, v2 = face[j], face[j+1]
                        p0, _ = to_screen_idx(v0)
                        p1, _ = to_screen_idx(v1)
                        p2, _ = to_screen_idx(v2)
                        if not (p0 and p1 and p2):
                            continue

                        if colorPerVertex and color and colorIndex:
                            # pega cor por vértice (em 0..1 → 0..255)
                            def col_of_vertex(ci):
                                # usa o índice correspondente; se faltar, usa emissive
                                if ci < len(colorIndex) and colorIndex[ci] >= 0:
                                    k = colorIndex[ci]*3
                                    rgb = [color[k]*255, color[k+1]*255, color[k+2]*255]
                                    return [int(max(0,min(255,round(v)))) for v in rgb]
                                return solid
                            c0 = col_of_vertex(v0)
                            c1 = col_of_vertex(v1)
                            c2 = col_of_vertex(v2)
                            fill_triangle_color(p0, p1, p2, c0, c1, c2)
                        else:
                            # cor sólida
                            GL._fill_triangle_simple(p0[0], p0[1], p1[0], p1[1], p2[0], p2[1], solid)
                face = []
            else:
                face.append(idx)

    @staticmethod
    def box(size, colors):
        """Função usada para renderizar Boxes."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Box
        # A função box é usada para desenhar paralelepípedos na cena. O Box é centrada no
        # (0, 0, 0) no sistema de coordenadas local e alinhado com os eixos de coordenadas
        # locais. O argumento size especifica as extensões da caixa ao longo dos eixos X, Y
        # e Z, respectivamente, e cada valor do tamanho deve ser maior que zero. Para desenha
        # essa caixa você vai provavelmente querer tesselar ela em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Box : size = {0}".format(size)) # imprime no terminal pontos
        print("Box : colors = {0}".format(colors)) # imprime no terminal as cores

        # Exemplo de desenho de um pixel branco na coordenada 10, 10
        gpu.GPU.draw_pixel([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

    @staticmethod
    def sphere(radius, colors):
        """Função usada para renderizar Esferas."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Sphere
        # A função sphere é usada para desenhar esferas na cena. O esfera é centrada no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da esfera que está sendo criada. Para desenha essa esfera você vai
        # precisar tesselar ela em triângulos, para isso encontre os vértices e defina
        # os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Sphere : radius = {0}".format(radius)) # imprime no terminal o raio da esfera
        print("Sphere : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def cone(bottomRadius, height, colors):
        """Função usada para renderizar Cones."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cone
        # A função cone é usada para desenhar cones na cena. O cone é centrado no
        # (0, 0, 0) no sistema de coordenadas local. O argumento bottomRadius especifica o
        # raio da base do cone e o argumento height especifica a altura do cone.
        # O cone é alinhado com o eixo Y local. O cone é fechado por padrão na base.
        # Para desenha esse cone você vai precisar tesselar ele em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Cone : bottomRadius = {0}".format(bottomRadius)) # imprime no terminal o raio da base do cone
        print("Cone : height = {0}".format(height)) # imprime no terminal a altura do cone
        print("Cone : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def cylinder(radius, height, colors):
        """Função usada para renderizar Cilindros."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cylinder
        # A função cylinder é usada para desenhar cilindros na cena. O cilindro é centrado no
        # (0, 0, 0) no sistema de coordenadas local. O argumento radius especifica o
        # raio da base do cilindro e o argumento height especifica a altura do cilindro.
        # O cilindro é alinhado com o eixo Y local. O cilindro é fechado por padrão em ambas as extremidades.
        # Para desenha esse cilindro você vai precisar tesselar ele em triângulos, para isso
        # encontre os vértices e defina os triângulos.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Cylinder : radius = {0}".format(radius)) # imprime no terminal o raio do cilindro
        print("Cylinder : height = {0}".format(height)) # imprime no terminal a altura do cilindro
        print("Cylinder : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def navigationInfo(headlight):
        """Características físicas do avatar do visualizador e do modelo de visualização."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/navigation.html#NavigationInfo
        # O campo do headlight especifica se um navegador deve acender um luz direcional que
        # sempre aponta na direção que o usuário está olhando. Definir este campo como TRUE
        # faz com que o visualizador forneça sempre uma luz do ponto de vista do usuário.
        # A luz headlight deve ser direcional, ter intensidade = 1, cor = (1 1 1),
        # ambientIntensity = 0,0 e direção = (0 0 −1).

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("NavigationInfo : headlight = {0}".format(headlight)) # imprime no terminal

    @staticmethod
    def directionalLight(ambientIntensity, color, intensity, direction):
        """Luz direcional ou paralela."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/lighting.html#DirectionalLight
        # Define uma fonte de luz direcional que ilumina ao longo de raios paralelos
        # em um determinado vetor tridimensional. Possui os campos básicos ambientIntensity,
        # cor, intensidade. O campo de direção especifica o vetor de direção da iluminação
        # que emana da fonte de luz no sistema de coordenadas local. A luz é emitida ao
        # longo de raios paralelos de uma distância infinita.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("DirectionalLight : ambientIntensity = {0}".format(ambientIntensity))
        print("DirectionalLight : color = {0}".format(color)) # imprime no terminal
        print("DirectionalLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("DirectionalLight : direction = {0}".format(direction)) # imprime no terminal

    @staticmethod
    def pointLight(ambientIntensity, color, intensity, location):
        """Luz pontual."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/lighting.html#PointLight
        # Fonte de luz pontual em um local 3D no sistema de coordenadas local. Uma fonte
        # de luz pontual emite luz igualmente em todas as direções; ou seja, é omnidirecional.
        # Possui os campos básicos ambientIntensity, cor, intensidade. Um nó PointLight ilumina
        # a geometria em um raio de sua localização. O campo do raio deve ser maior ou igual a
        # zero. A iluminação do nó PointLight diminui com a distância especificada.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("PointLight : ambientIntensity = {0}".format(ambientIntensity))
        print("PointLight : color = {0}".format(color)) # imprime no terminal
        print("PointLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("PointLight : location = {0}".format(location)) # imprime no terminal

    @staticmethod
    def fog(visibilityRange, color):
        """Névoa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/environmentalEffects.html#Fog
        # O nó Fog fornece uma maneira de simular efeitos atmosféricos combinando objetos
        # com a cor especificada pelo campo de cores com base nas distâncias dos
        # vários objetos ao visualizador. A visibilidadeRange especifica a distância no
        # sistema de coordenadas local na qual os objetos são totalmente obscurecidos
        # pela névoa. Os objetos localizados fora de visibilityRange do visualizador são
        # desenhados com uma cor de cor constante. Objetos muito próximos do visualizador
        # são muito pouco misturados com a cor do nevoeiro.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("Fog : color = {0}".format(color)) # imprime no terminal
        print("Fog : visibilityRange = {0}".format(visibilityRange))

    @staticmethod
    def timeSensor(cycleInterval, loop):
        """Gera eventos conforme o tempo passa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/time.html#TimeSensor
        # Os nós TimeSensor podem ser usados para muitas finalidades, incluindo:
        # Condução de simulações e animações contínuas; Controlar atividades periódicas;
        # iniciar eventos de ocorrência única, como um despertador;
        # Se, no final de um ciclo, o valor do loop for FALSE, a execução é encerrada.
        # Por outro lado, se o loop for TRUE no final de um ciclo, um nó dependente do
        # tempo continua a execução no próximo ciclo. O ciclo de um nó TimeSensor dura
        # cycleInterval segundos. O valor de cycleInterval deve ser maior que zero.

        # Deve retornar a fração de tempo passada em fraction_changed

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("TimeSensor : cycleInterval = {0}".format(cycleInterval)) # imprime no terminal
        print("TimeSensor : loop = {0}".format(loop))

        # Esse método já está implementado para os alunos como exemplo
        epoch = time.time()  # time in seconds since the epoch as a floating point number.
        fraction_changed = (epoch % cycleInterval) / cycleInterval

        return fraction_changed

    @staticmethod
    def splinePositionInterpolator(set_fraction, key, keyValue, closed):
        """Interpola não linearmente entre uma lista de vetores 3D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#SplinePositionInterpolator
        # Interpola não linearmente entre uma lista de vetores 3D. O campo keyValue possui
        # uma lista com os valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantos vetores 3D quanto os
        # quadros-chave no key. O campo closed especifica se o interpolador deve tratar a malha
        # como fechada, com uma transições da última chave para a primeira chave. Se os keyValues
        # na primeira e na última chave não forem idênticos, o campo closed será ignorado.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("SplinePositionInterpolator : set_fraction = {0}".format(set_fraction))
        print("SplinePositionInterpolator : key = {0}".format(key)) # imprime no terminal
        print("SplinePositionInterpolator : keyValue = {0}".format(keyValue))
        print("SplinePositionInterpolator : closed = {0}".format(closed))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0.0, 0.0, 0.0]
        
        return value_changed

    @staticmethod
    def orientationInterpolator(set_fraction, key, keyValue):
        """Interpola entre uma lista de valores de rotação especificos."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#OrientationInterpolator
        # Interpola rotações são absolutas no espaço do objeto e, portanto, não são cumulativas.
        # Uma orientação representa a posição final de um objeto após a aplicação de uma rotação.
        # Um OrientationInterpolator interpola entre duas orientações calculando o caminho mais
        # curto na esfera unitária entre as duas orientações. A interpolação é linear em
        # comprimento de arco ao longo deste caminho. Os resultados são indefinidos se as duas
        # orientações forem diagonalmente opostas. O campo keyValue possui uma lista com os
        # valores a serem interpolados, key possui uma lista respectiva de chaves
        # dos valores em keyValue, a fração a ser interpolada vem de set_fraction que varia de
        # zeroa a um. O campo keyValue deve conter exatamente tantas rotações 3D quanto os
        # quadros-chave no key.

        # O print abaixo é só para vocês verificarem o funcionamento, DEVE SER REMOVIDO.
        print("OrientationInterpolator : set_fraction = {0}".format(set_fraction))
        print("OrientationInterpolator : key = {0}".format(key)) # imprime no terminal
        print("OrientationInterpolator : keyValue = {0}".format(keyValue))

        # Abaixo está só um exemplo de como os dados podem ser calculados e transferidos
        value_changed = [0, 0, 1, 0]

        return value_changed

    # Para o futuro (Não para versão atual do projeto.)
    def vertex_shader(self, shader):
        """Para no futuro implementar um vertex shader."""

    def fragment_shader(self, shader):
        """Para no futuro implementar um fragment shader."""

    # Funções auxiliares simples para rasterização
    @staticmethod
    def _draw_line_simple(x0, y0, x1, y1, color):
        """Desenha uma linha simples entre dois pontos."""
        # Algoritmo simples de linha
        dx = abs(x1 - x0)
        dy = abs(y1 - y0)
        
        if dx > dy:
            # Linha mais horizontal
            if x0 > x1:
                x0, x1 = x1, x0
                y0, y1 = y1, y0
            
            for x in range(x0, x1 + 1):
                y = int(y0 + (y1 - y0) * (x - x0) / (x1 - x0))
                if 0 <= x < GL.width and 0 <= y < GL.height:
                    gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, color)
        else:
            # Linha mais vertical
            if y0 > y1:
                x0, x1 = x1, x0
                y0, y1 = y1, y0
            
            for y in range(y0, y1 + 1):
                x = int(x0 + (x1 - x0) * (y - y0) / (y1 - y0))
                if 0 <= x < GL.width and 0 <= y < GL.height:
                    gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, color)

    @staticmethod
    def _fill_triangle_simple(x1, y1, x2, y2, x3, y3, color):
        """Preenche um triângulo de forma simples."""
        # Ordenar vértices por y (do menor para o maior)
        vertices = [(x1, y1), (x2, y2), (x3, y3)]
        vertices.sort(key=lambda v: v[1])
        
        x1, y1 = vertices[0]
        x2, y2 = vertices[1]
        x3, y3 = vertices[2]
        
        # Preencher triângulo usando algoritmo de scanline simples
        for y in range(int(y1), int(y3) + 1):
            if 0 <= y < GL.height:
                # Calcular x inicial e final para esta linha y
                if y <= y2:
                    # Primeira parte do triângulo (y1 até y2)
                    if y2 != y1:
                        x_start = int(x1 + (x2 - x1) * (y - y1) / (y2 - y1))
                        x_end = int(x1 + (x3 - x1) * (y - y1) / (y3 - y1))
                    else:
                        x_start = x1
                        x_end = x1
                else:
                    # Segunda parte do triângulo (y2 até y3)
                    if y3 != y2:
                        x_start = int(x2 + (x3 - x2) * (y - y2) / (y3 - y2))
                        x_end = int(x1 + (x3 - x1) * (y - y1) / (y3 - y1))
                    else:
                        x_start = x2
                        x_end = x2
                
                # Ordenar x_start e x_end
                if x_start > x_end:
                    x_start, x_end = x_end, x_start
                
                # Desenhar linha horizontal preenchendo o triângulo
                for x in range(x_start, x_end + 1):
                    if 0 <= x < GL.width:
                        gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, color)
