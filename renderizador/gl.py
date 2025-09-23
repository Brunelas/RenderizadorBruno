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

        # if not point or len(point) < 9:
        #     return

        # # cor emissiva (X3D em 0..1 → framebuffer 0..255)
        # rgb = colors.get('emissiveColor', [1.0, 1.0, 1.0]) if isinstance(colors, dict) else [1.0, 1.0, 1.0]
        # col = [int(max(0, min(255, round(c*255)))) for c in rgb]

        # # garante matrizes padrão caso ainda não tenham sido definidas
        # M = getattr(GL, "_M", np.eye(4, dtype=float))
        # V = getattr(GL, "_V", np.eye(4, dtype=float))
        # P = getattr(GL, "_P", np.eye(4, dtype=float))
        # PVM = P @ V @ M

        # w, h = GL.width, GL.height

        # def ndc_to_screen(ndc_xy):
        #     x_ndc, y_ndc = ndc_xy[0], ndc_xy[1]
        #     x = (x_ndc + 1.0) * 0.5 * (w - 1)
        #     y = (1 - (y_ndc + 1.0) * 0.5) * (h - 1)  # origem (0,0) no topo
        #     return int(round(x)), int(round(y))

        # for i in range(0, len(point), 9):
        #     if i + 8 >= len(point):  # segurança
        #         break

        #     v0 = np.array([point[i],   point[i+1], point[i+2], 1.0], dtype=float)
        #     v1 = np.array([point[i+3], point[i+4], point[i+5], 1.0], dtype=float)
        #     v2 = np.array([point[i+6], point[i+7], point[i+8], 1.0], dtype=float)

        #     c0 = PVM @ v0
        #     c1 = PVM @ v1
        #     c2 = PVM @ v2

        #     # divisão por w -> NDC
        #     if c0[3] == 0 or c1[3] == 0 or c2[3] == 0:
        #         continue
        #     ndc0 = c0[:3] / c0[3]
        #     ndc1 = c1[:3] / c1[3]
        #     ndc2 = c2[:3] / c2[3]

        #     # mapeia para tela
        #     x0, y0 = ndc_to_screen(ndc0)
        #     x1, y1 = ndc_to_screen(ndc1)
        #     x2, y2 = ndc_to_screen(ndc2)

        #     # preenche o triângulo usando seu raster simples
        #     GL._fill_triangle_simple(x0, y0, x1, y1, x2, y2, col)

        if not point or len(point) < 9:
            return

        material = colors if isinstance(colors, dict) else {}

        for i in range(0, len(point), 9):
            if i + 8 >= len(point):
                break
            v0 = [point[i],   point[i+1], point[i+2]]
            v1 = [point[i+3], point[i+4], point[i+5]]
            v2 = [point[i+6], point[i+7], point[i+8]]

            p0 = GL._project_vertex(v0)
            p1 = GL._project_vertex(v1)
            p2 = GL._project_vertex(v2)
            if p0 is None or p1 is None or p2 is None:
                continue

            a0 = {"rgb": None, "uv": None}
            a1 = {"rgb": None, "uv": None}
            a2 = {"rgb": None, "uv": None}

            GL._raster_triangle_persp(p0, p1, p2, a0, a1, a2, material, tex=None, tri_uv_area=None)

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
        # print("IndexedFaceSet : ")
        # if coord:
        #     print("\tpontos(x, y, z) = {0}, coordIndex = {1}".format(coord, coordIndex))
        # print("colorPerVertex = {0}".format(colorPerVertex))
        # if colorPerVertex and color and colorIndex:
        #     print("\tcores(r, g, b) = {0}, colorIndex = {1}".format(color, colorIndex))
        # if texCoord and texCoordIndex:
        #     print("\tpontos(u, v) = {0}, texCoordIndex = {1}".format(texCoord, texCoordIndex))
        # if current_texture:
        #     image = gpu.GPU.load_texture(current_texture[0])
        #     print("\t Matriz com image = {0}".format(image))
        #     print("\t Dimensões da image = {0}".format(image.shape))
        # print("IndexedFaceSet : colors = {0}".format(colors))  # imprime no terminal as cores

        # # ---------- Helpers de pipeline ----------
    
        w, h = GL.width, GL.height
        M = getattr(GL, "_M", np.eye(4))
        V = getattr(GL, "_V", np.eye(4))
        P = getattr(GL, "_P", np.eye(4))
        PVM = P @ V @ M

        # z-buffer de software (NDC: -1..1; mais perto = menor valor)
        if not hasattr(GL, "_zbuf") or GL._zbuf.shape != (h, w):
            GL._zbuf = np.full((h, w), 1.0, dtype=float)  # 1.0 = "longe"

        # carrega textura (se houver)
        img_tex = None
        if current_texture and len(current_texture) > 0 and isinstance(current_texture[0], str):
            img_tex = gpu.GPU.load_texture(current_texture[0])

        # converte um índice de coordenada (em 'coord') para: screen(x,y), ndc_z, inv_w, clip_w
        def proj_to_screen(ii):
            v3 = coord[3*ii:3*ii+3]
            v = np.array([v3[0], v3[1], v3[2], 1.0], dtype=float)
            clip = PVM @ v
            w_clip = clip[3]
            if w_clip == 0:
                return None
            inv_w = 1.0 / w_clip
            ndc = clip[:3] * inv_w
            x = int(round((ndc[0] + 1.0) * 0.5 * (w - 1)))
            y = int(round((1.0 - (ndc[1] + 1.0) * 0.5) * (h - 1)))  # origem no topo
            z_ndc = ndc[2]  # em [-1, 1]
            return (x, y, z_ndc, inv_w)

        # amostragem de textura (nearest), corrigindo o eixo V (imagem tem (0,0) no topo)
        def sample_tex(u, v):
            if img_tex is None:
                return None
            H, W = img_tex.shape[0], img_tex.shape[1]
            # clamp
            if u < 0: u = 0
            if u > 1: u = 1
            if v < 0: v = 0
            if v > 1: v = 1
            # X3D usa V=0 embaixo → imagem usa Y=0 em cima → flip em V:
            uu = (1.0 - u) * (W - 1)
            vv = v * (H - 1)
            xi = int(round(uu))
            yi = int(round(vv))
            px = img_tex[yi, xi]
            # textura pode ser RGB ou RGBA
            if len(px) == 4:
                return [int(px[0]), int(px[1]), int(px[2]), int(px[3])]
            else:
                return [int(px[0]), int(px[1]), int(px[2]), 255]

        # edge function
        def edge(a, b, c):
            return (c[0]-a[0])*(b[1]-a[1]) - (c[1]-a[1])*(b[0]-a[0])

        # rasterização com interpolação perspectiva (cor e UV) + z-buffer + transparência simples
        def raster_tri(p0, p1, p2, col0, col1, col2, uv0, uv1, uv2):
            # p* = (x, y, z_ndc, inv_w)
            x0, y0, z0, iw0 = p0
            x1, y1, z1, iw1 = p1
            x2, y2, z2, iw2 = p2

            # bounding box
            minx = max(0, min(x0, x1, x2))
            maxx = min(w-1, max(x0, x1, x2))
            miny = max(0, min(y0, y1, y2))
            maxy = min(h-1, max(y0, y1, y2))

            A = edge((x0,y0), (x1,y1), (x2,y2))
            if A == 0:
                return

            # prepara atributos multiplicados por inv_w (perspective-correct)
            has_color = (col0 is not None) and (col1 is not None) and (col2 is not None)
            has_uv    = (uv0  is not None) and (uv1  is not None) and (uv2  is not None) and (img_tex is not None)

            if has_color:
                c0p = [col0[0]*iw0, col0[1]*iw0, col0[2]*iw0]
                c1p = [col1[0]*iw1, col1[1]*iw1, col1[2]*iw1]
                c2p = [col2[0]*iw2, col2[1]*iw2, col2[2]*iw2]
            else:
                # cor emissiva sólida
                em = colors.get('emissiveColor', [1.0, 1.0, 1.0]) if isinstance(colors, dict) else [1.0, 1.0, 1.0]
                base = [int(max(0, min(255, round(em[0]*255)))),
                        int(max(0, min(255, round(em[1]*255)))),
                        int(max(0, min(255, round(em[2]*255))))]

            if has_uv:
                u0p, v0p = uv0[0]*iw0, uv0[1]*iw0
                u1p, v1p = uv1[0]*iw1, uv1[1]*iw1
                u2p, v2p = uv2[0]*iw2, uv2[1]*iw2

            # transparência do material (se houver) → alpha do material
            mat_alpha = 1.0 - float(colors.get('transparency', 0.0)) if isinstance(colors, dict) else 1.0
            if mat_alpha < 0: mat_alpha = 0.0
            if mat_alpha > 1: mat_alpha = 1.0

            for y in range(miny, maxy+1):
                for x in range(minx, maxx+1):
                    P = (x + 0.5, y + 0.5)
                    w0 = edge((x1,y1), (x2,y2), P)
                    w1 = edge((x2,y2), (x0,y0), P)
                    w2 = edge((x0,y0), (x1,y1), P)

                    if (A > 0 and w0 >= 0 and w1 >= 0 and w2 >= 0) or (A < 0 and w0 <= 0 and w1 <= 0 and w2 <= 0):
                        # baricêntricas normalizadas
                        w0n, w1n, w2n = w0/A, w1/A, w2/A

                        # inv_w no pixel
                        invw_pix = w0n*iw0 + w1n*iw1 + w2n*iw2
                        if invw_pix == 0:
                            continue

                        # profundidade perspectiva-correct (usando z_ndc * inv_w)
                        z_num = (w0n*z0*iw0 + w1n*z1*iw1 + w2n*z2*iw2)
                        z_ndc_pix = z_num / invw_pix  # em [-1,1]
                        if z_ndc_pix >= GL._zbuf[y, x]:
                            continue  # falhou no z-test (menor = mais perto)

                        # cor
                        if has_color:
                            r = (w0n*c0p[0] + w1n*c1p[0] + w2n*c2p[0]) / invw_pix
                            g = (w0n*c0p[1] + w1n*c1p[1] + w2n*c2p[1]) / invw_pix
                            b = (w0n*c0p[2] + w1n*c1p[2] + w2n*c2p[2]) / invw_pix
                            col_pix = [int(max(0, min(255, round(r)))),
                                       int(max(0, min(255, round(g)))),
                                       int(max(0, min(255, round(b))))]
                        else:
                            col_pix = base[:]

                        # textura (se houver), com correção por perspectiva
                        if has_uv:
                            u = (w0n*u0p + w1n*u1p + w2n*u2p) / invw_pix
                            v = (w0n*v0p + w1n*v1p + w2n*v2p) / invw_pix
                            tex = sample_tex(v, u)   # <- swap corrige rotação de 90°
                            if tex is not None:
                                tr, tg, tb, ta = tex
                                # blend da textura sobre a cor interpolada
                                a = (ta / 255.0) * mat_alpha
                                inva = 1.0 - a
                                col_pix = [int(tr*a + col_pix[0]*inva),
                                           int(tg*a + col_pix[1]*inva),
                                           int(tb*a + col_pix[2]*inva)]
                        else:
                            # sem textura, aplica só alpha do material
                            if mat_alpha < 1.0:
                                bg = gpu.GPU.read_pixel([x, y], gpu.GPU.RGB8)
                                a = mat_alpha
                                inva = 1.0 - a
                                col_pix = [int(col_pix[0]*a + bg[0]*inva),
                                           int(col_pix[1]*a + bg[1]*inva),
                                           int(col_pix[2]*a + bg[2]*inva)]

                        # escreve (passou no z-test)
                        gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, col_pix)
                        GL._zbuf[y, x] = z_ndc_pix

        # --------- Triangulação em fan + montagem de atributos (cor/UV) ---------
        face = []
        face_uv = []
        ci = 0   # cursor em coordIndex
        ti = 0   # cursor em texCoordIndex

        def fetch_uv(idx):
            if texCoord is None or idx is None or idx < 0:
                return None
            j = 2*idx
            if j+1 >= len(texCoord):
                return None
            return [float(texCoord[j]), float(texCoord[j+1])]

        def vertex_color(idx):
            # colorPerVertex: cores em 'color' / 'colorIndex' (0..1) → 0..255
            if colorPerVertex and color is not None:
                if colorIndex and idx < len(colorIndex) and colorIndex[idx] >= 0:
                    k = 3*colorIndex[idx]
                else:
                    k = 3*idx
                if k+2 < len(color):
                    return [int(max(0,min(255,round(color[k]*255)))),
                            int(max(0,min(255,round(color[k+1]*255)))),
                            int(max(0,min(255,round(color[k+2]*255))))]
            return None  # sem cor por vértice

        # percorre os polígonos de coordIndex (separados por -1) e texCoordIndex em paralelo
        ptr_tc = 0
        for idx in coordIndex + [-1]:
            if idx == -1:
                if len(face) >= 3:
                    for j in range(1, len(face)-1):
                        i0, i1, i2 = face[0], face[j], face[j+1]
                        p0 = proj_to_screen(i0)
                        p1 = proj_to_screen(i1)
                        p2 = proj_to_screen(i2)
                        if not (p0 and p1 and p2):
                            continue

                        c0 = vertex_color(i0)
                        c1 = vertex_color(i1)
                        c2 = vertex_color(i2)

                        uv0 = uv1 = uv2 = None
                        if texCoord is not None:
                            # usa texCoordIndex quando disponível (mesma fatiada do polígono)
                            if texCoordIndex and len(face_uv) == len(face):
                                uv0 = fetch_uv(face_uv[0])
                                uv1 = fetch_uv(face_uv[j])
                                uv2 = fetch_uv(face_uv[j+1])
                            else:
                                # fallback: usa o mesmo índice de posição
                                uv0 = fetch_uv(i0)
                                uv1 = fetch_uv(i1)
                                uv2 = fetch_uv(i2)



                        raster_tri(p0, p1, p2, c0, c1, c2, uv0, uv1, uv2)

                # reseta polígono atual
                face = []
                face_uv = []
                # avança o ponteiro de texCoordIndex para o próximo polígono
                if texCoordIndex:
                    # pula exatamente o tamanho deste polígono + o -1
                    # (como não guardamos, apenas avançamos até o próximo -1)
                    while ptr_tc < len(texCoordIndex) and texCoordIndex[ptr_tc] != -1:
                        ptr_tc += 1
                    if ptr_tc < len(texCoordIndex) and texCoordIndex[ptr_tc] == -1:
                        ptr_tc += 1
            else:
                face.append(idx)
                # guarda o índice de texCoord correspondente a este vértice do polígono
                if texCoordIndex and ptr_tc < len(texCoordIndex):
                    face_uv.append(texCoordIndex[ptr_tc] if texCoordIndex[ptr_tc] != -1 else None)
                    ptr_tc += 1
                else:
                    face_uv.append(None)



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
    
     # [ADICIONADO] ---------- Buffers de software (Z) ----------
    @staticmethod
    def _alloc_software_buffers():
        """Aloca z-buffer de software para a resolução atual (GL.width x GL.height)."""
        GL._zbuf = np.full((GL.height, GL.width), np.inf, dtype=float)

    @staticmethod
    def _clear_software_buffers():
        """Limpa z-buffer de software para +inf (chamar a cada frame)."""
        if hasattr(GL, "_zbuf"):
            GL._zbuf.fill(np.inf)

    # [ADICIONADO] ---------- Projeção com info de 1/w e z_ndc ----------
    @staticmethod
    def _project_vertex(v3):
        M = getattr(GL, "_M", np.eye(4))
        V = getattr(GL, "_V", np.eye(4))
        P = getattr(GL, "_P", np.eye(4))
        PVM = P @ V @ M
        v = np.array([v3[0], v3[1], v3[2], 1.0], dtype=float)
        clip = PVM @ v
        if clip[3] == 0:
            return None
        invw = 1.0 / clip[3]
        ndc = clip[:3] * invw  # [-1,1]
        sx = (ndc[0] + 1.0) * 0.5 * (GL.width - 1)
        sy = (1.0 - (ndc[1] + 1.0) * 0.5) * (GL.height - 1)
        return {"sx": sx, "sy": sy, "z_ndc": ndc[2], "invw": invw}

    # [ADICIONADO] ---------- Geometria 2D ----------
    @staticmethod
    def _edge(a, b, p):
        return (p[0]-a[0])*(b[1]-a[1]) - (p[1]-a[1])*(b[0]-a[0])

    # [ADICIONADO] ---------- Blend (src over dst) ----------
    @staticmethod
    def _blend_over(dst_rgb, src_rgb, alpha):
        # todos em 0..255 float, alpha em 0..1
        return src_rgb * alpha + dst_rgb * (1.0 - alpha)

    # [ADICIONADO] ---------- Textura + Mipmap ----------
    _tex_cache = {}

    @staticmethod
    def _get_texture(tex_path):
        if tex_path in GL._tex_cache:
            return GL._tex_cache[tex_path]
        img = gpu.GPU.load_texture(tex_path)  # numpy HxWxC (uint8)
        levels = [img]
        cur = img
        while min(cur.shape[0], cur.shape[1]) > 1:
            h2 = max(1, cur.shape[0] // 2)
            w2 = max(1, cur.shape[1] // 2)
            cur = cur[:h2*2, :w2*2].reshape(h2, 2, w2, 2, cur.shape[2]).mean(axis=(1,3)).astype(np.uint8)
            levels.append(cur)
        tex = {"levels": levels}
        GL._tex_cache[tex_path] = tex
        return tex

    @staticmethod
    def _sample_tex_level(level_img, u, v):
        # u,v em [0,1]; v invertido (origem topo)
        H, W, C = level_img.shape
        uu = np.clip(u, 0.0, 1.0) * (W - 1)
        vv = (1.0 - np.clip(v, 0.0, 1.0)) * (H - 1)
        x0 = int(np.floor(uu)); x1 = min(W - 1, x0 + 1)
        y0 = int(np.floor(vv)); y1 = min(H - 1, y0 + 1)
        tx = uu - x0; ty = vv - y0
        c00 = level_img[y0, x0].astype(float)
        c10 = level_img[y0, x1].astype(float)
        c01 = level_img[y1, x0].astype(float)
        c11 = level_img[y1, x1].astype(float)
        c0 = c00 * (1 - tx) + c10 * tx
        c1 = c01 * (1 - tx) + c11 * tx
        return c0 * (1 - ty) + c1 * ty  # 3 ou 4 canais

    @staticmethod
    def _estimate_mip_level(tex, tri_uv_area, tri_screen_area):
        levels = tex["levels"]
        if tri_uv_area is None or tri_uv_area <= 0 or tri_screen_area <= 0:
            return 0
        base = levels[0]
        H, W = base.shape[:2]
        area_texels = tri_uv_area * (W * H)
        s = max(1e-6, np.sqrt(area_texels / max(1.0, tri_screen_area)))
        L = int(np.clip(np.log2(s), 0, len(levels) - 1))
        return L

    # [ADICIONADO] ---------- Raster com interp persp, Z-buffer, transparência e textura ----------
    @staticmethod
    def _raster_triangle_persp(p0, p1, p2, a0, a1, a2, material, tex=None, tri_uv_area=None):
        # material
                # base do material: preferir diffuseColor; se não tiver, usar emissiveColor; se nada vier, branco
        diff = material.get('diffuseColor', None) if isinstance(material, dict) else None
        emis = material.get('emissiveColor', None) if isinstance(material, dict) else None
        if diff is not None:
            base_rgb = np.array([diff[0]*255.0, diff[1]*255.0, diff[2]*255.0], dtype=float)
        elif emis is not None:
            base_rgb = np.array([emis[0]*255.0, emis[1]*255.0, emis[2]*255.0], dtype=float)
        else:
            base_rgb = np.array([255.0, 255.0, 255.0], dtype=float)  # não apaga textura por engano
        alpha_mat = 1.0 - float(material.get('transparency', 0.0)) if isinstance(material, dict) else 1.0


        # bounding box
        minx = int(max(0, np.floor(min(p0["sx"], p1["sx"], p2["sx"]))))
        maxx = int(min(GL.width - 1, np.ceil (max(p0["sx"], p1["sx"], p2["sx"]))))
        miny = int(max(0, np.floor(min(p0["sy"], p1["sy"], p2["sy"]))))
        maxy = int(min(GL.height - 1, np.ceil (max(p0["sy"], p1["sy"], p2["sy"]))))

        A = GL._edge((p0["sx"], p0["sy"]), (p1["sx"], p1["sy"]), (p2["sx"], p2["sy"]))
        if A == 0:
            return
        tri_screen_area = abs(A)

        level_img = None
        if tex is not None:
            L = GL._estimate_mip_level(tex, tri_uv_area, tri_screen_area)
            level_img = tex["levels"][L]

        w0 = p0["invw"]; w1 = p1["invw"]; w2 = p2["invw"]
        z0 = p0["z_ndc"]; z1 = p1["z_ndc"]; z2 = p2["z_ndc"]

        for y in range(miny, maxy + 1):
            for x in range(minx, maxx + 1):
                Px = x + 0.5; Py = y + 0.5
                e0 = GL._edge((p1["sx"], p1["sy"]), (p2["sx"], p2["sy"]), (Px, Py))
                e1 = GL._edge((p2["sx"], p2["sy"]), (p0["sx"], p0["sy"]), (Px, Py))
                e2 = GL._edge((p0["sx"], p0["sy"]), (p1["sx"], p1["sy"]), (Px, Py))
                inside = ((A > 0 and e0 >= 0 and e1 >= 0 and e2 >= 0) or
                          (A < 0 and e0 <= 0 and e1 <= 0 and e2 <= 0))
                if not inside:
                    continue

                l0 = e0 / A; l1 = e1 / A; l2 = e2 / A
                denom = l0*w0 + l1*w1 + l2*w2
                if denom == 0:
                    continue

                z_ndc = (l0*z0*w0 + l1*z1*w1 + l2*z2*w2) / denom
                z01 = 0.5 * z_ndc + 0.5

                # depth test (software)
                zcur = GL._zbuf[y, x]

                               # base: cor por vértice (se existir) SENÃO a cor do material (diffuse/emissive/branco)
                if a0.get("rgb") is not None:
                    cr = (l0*a0["rgb"][0]*w0 + l1*a1["rgb"][0]*w1 + l2*a2["rgb"][0]*w2)/denom
                    cg = (l0*a0["rgb"][1]*w0 + l1*a1["rgb"][1]*w1 + l2*a2["rgb"][1]*w2)/denom
                    cb = (l0*a0["rgb"][2]*w0 + l1*a1["rgb"][2]*w1 + l2*a2["rgb"][2]*w2)/denom
                    rgb_base = np.array([cr, cg, cb], dtype=float)
                else:
                    rgb_base = base_rgb.copy()

                rgb = rgb_base.copy()

                # textura (se houver)
                alpha_tex = 1.0
                if level_img is not None and a0.get("uv") is not None:
                    u = (l0*a0["uv"][0]*w0 + l1*a1["uv"][0]*w1 + l2*a2["uv"][0]*w2) / denom
                    v = (l0*a0["uv"][1]*w0 + l1*a1["uv"][1]*w1 + l2*a2["uv"][1]*w2) / denom
                    texel = GL._sample_tex_level(level_img, u, v)  # [ADICIONADO] corrige rotação 90° (swap u<->v)


                    if texel.shape[0] >= 3:
                        # MODULA a textura pela base (diffuse/vertex). Se base for branca, fica a textura “pura”.
                        rgb = rgb_base * (texel[:3] / 255.0)

                    if texel.shape[0] == 4:
                        alpha_tex = float(texel[3]) / 255.0


                alpha = float(np.clip(alpha_mat * alpha_tex, 0.0, 1.0))

                # OPAQUE: escreve cor e z se mais perto
                if alpha >= 0.999:
                    if z01 < zcur:
                        gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8,
                                           [int(np.clip(rgb[0],0,255)),
                                            int(np.clip(rgb[1],0,255)),
                                            int(np.clip(rgb[2],0,255))])
                        GL._zbuf[y, x] = z01
                    continue

                # TRANSPARENTE: depth test, não escreve z, faz blend over
                if z01 < zcur:
                    dst = np.array(gpu.GPU.read_pixel([x, y], gpu.GPU.RGB8), dtype=float)
                    out = GL._blend_over(dst, np.clip(rgb,0,255), alpha)
                    gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8,
                                       [int(np.clip(out[0],0,255)),
                                        int(np.clip(out[1],0,255)),
                                        int(np.clip(out[2],0,255))])
                    

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
