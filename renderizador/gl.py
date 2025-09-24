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
    # [ADICIONADO iluminação] estado global de luzes
    _lights = []                 # cada item: {"dir": [dx,dy,dz], "color":[r,g,b], "intensity": f, "ambient": f}
    _headlight_on = True         # X3D default = TRUE
    # estado de iluminação
    _directional_light = None   # dict com a luz direcional atual

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
        """Função usada para renderizar TriangleSet (sem iluminação: cor chapada)."""
        if not point or len(point) < 9:
            return

        material = colors if isinstance(colors, dict) else {}

        for i in range(0, len(point), 9):
            if i + 8 >= len(point):
                break

            v0 = [point[i],   point[i+1], point[i+2]]
            v1 = [point[i+3], point[i+4], point[i+5]]
            v2 = [point[i+6], point[i+7], point[i+8]]

            # projeta para tela (P*V*M já está embutido em _project_vertex)
            p0 = GL._project_vertex(v0)
            p1 = GL._project_vertex(v1)
            p2 = GL._project_vertex(v2)
            if p0 is None or p1 is None or p2 is None:
                continue

            # NÃO passar normal -> raster não aplica iluminação
            a0 = {"rgb": None, "uv": None, "n": None}
            a1 = {"rgb": None, "uv": None, "n": None}
            a2 = {"rgb": None, "uv": None, "n": None}

            GL._raster_triangle_persp(p0, p1, p2, a0, a1, a2, material, tex=None, tri_uv_area=None)


    @staticmethod
    def viewpoint(position, orientation, fieldOfView):
        """Função usada para renderizar (na verdade coletar os dados) de Viewpoint."""
        # Na função de viewpoint você receberá a posição, orientação e campo de visão da
        # câmera virtual.

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

        setattr(GL, "_V", V)   # (mantém só uma vez)
        setattr(GL, "_P", P)

        # headlight básico (só se ninguém criou outra DirectionalLight)
        if getattr(GL, "_headlight_on", True) and getattr(GL, "_directional_light", None) is None:
            GL._directional_light = {
                "ambientIntensity": 0.15,
                "intensity": 1.1,
                "color": [1.0, 1.0, 1.0],
                "direction_world": np.array([0.0, 0.0, -1.0], dtype=float),
            }

        # [ADICIONADO] guarda posição e direção da câmera para iluminação / headlight
        setattr(GL, "_cam_pos", np.array([px, py, pz], dtype=float))
        # direção -Z no espaço da câmera, rotacionada para mundo (headlight segue olhar)
        R_cam = R_axis_angle(ox, oy, oz, ang)
        fwd = (R_cam @ np.array([0.0, 0.0, -1.0, 0.0]))[:3]
        setattr(GL, "_cam_fwd", fwd / (np.linalg.norm(fwd) or 1.0))

    


    @staticmethod
    def transform_in(translation, scale, rotation):
        """Função usada para renderizar (na verdade coletar os dados) de Transform."""
        # A função transform_in será chamada quando se entrar em um nó X3D do tipo Transform

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
        # restaura a matriz do topo da pilha
        if hasattr(GL, "_M_stack") and len(GL._M_stack) > 0:
            GL._M = GL._M_stack.pop()
        else:
            # se algo sair do esperado, volta para identidade
            GL._M = np.eye(4, dtype=float)

    @staticmethod
    def triangleStripSet(point, stripCount, colors):
        """Desenha TriangleStripSet com sombreamento suave (normais por vértice)."""
        if not point or not stripCount:
            return

        # --- dados básicos ---
        P = np.array(point, dtype=float).reshape((-1, 3))  # Nx3 (world)
        nverts = P.shape[0]

        # mundo -> câmera (para computar normal no espaço de visão)
        V = getattr(GL, "_V", np.eye(4))
        M = getattr(GL, "_M", np.eye(4))
        MV = V @ M

        # posições no espaço da câmera (uma vez só)
        Pw = np.hstack([P, np.ones((nverts, 1), dtype=float)])         # Nx4
        Pe = (MV @ Pw.T).T[:, :3]                                      # Nx3 (eye space)

        # --- acumula normais por vértice ---
        N = np.zeros_like(Pe)  # Nx3
        base = 0
        for count in stripCount:
            idx = list(range(base, base + count))
            for i in range(count - 2):
                # alterna a orientação típica de strip
                if i % 2 == 0:
                    i0, i1, i2 = idx[i], idx[i+1], idx[i+2]
                else:
                    i0, i1, i2 = idx[i+1], idx[i], idx[i+2]

                e0, e1, e2 = Pe[i0], Pe[i1], Pe[i2]
                faceN = np.cross(e1 - e0, e2 - e0)
                lenN = np.linalg.norm(faceN)
                if lenN < 1e-12:
                    continue
                faceN /= lenN
                N[i0] += faceN; N[i1] += faceN; N[i2] += faceN
            base += count

        # normaliza por vértice
        lens = np.linalg.norm(N, axis=1)
        lens[lens == 0] = 1.0
        N /= lens[:, None]

        # --- rasteriza usando normais por vértice ---
        base = 0
        for count in stripCount:
            idx = list(range(base, base + count))
            for i in range(count - 2):
                if i % 2 == 0:
                    i0, i1, i2 = idx[i], idx[i+1], idx[i+2]
                else:
                    i0, i1, i2 = idx[i+1], idx[i], idx[i+2]

                v0 = P[i0].tolist(); v1 = P[i1].tolist(); v2 = P[i2].tolist()

                # projeção (P*V*M) -> já devolve sx, sy, invw, z_ndc etc.
                p0 = GL._project_vertex(v0)
                p1 = GL._project_vertex(v1)
                p2 = GL._project_vertex(v2)
                if p0 is None or p1 is None or p2 is None:
                    continue

                # atributos por-vértice: normal SUAVE (no espaço da câmera)
                a0 = {"rgb": None, "uv": None, "n": N[i0].tolist()}
                a1 = {"rgb": None, "uv": None, "n": N[i1].tolist()}
                a2 = {"rgb": None, "uv": None, "n": N[i2].tolist()}

                GL._raster_triangle_persp(p0, p1, p2, a0, a1, a2, colors, tex=None, tri_uv_area=None)

            base += count



    @staticmethod
    def indexedTriangleStripSet(point, index, colors):
        """Renderiza IndexedTriangleStripSet com projeção, z-buffer e iluminação."""
        if not point or not index:
            return

        M = getattr(GL, "_M", np.eye(4))
        V = getattr(GL, "_V", np.eye(4))
        MV = V @ M

        def vtx(i):
            j = 3*i
            return [float(point[j]), float(point[j+1]), float(point[j+2])]

        def to_eye(v):
            P = np.array([v[0], v[1], v[2], 1.0], dtype=float)
            e = MV @ P
            return e[:3]

        strip = []
        for idx in index + [-1]:
            if idx == -1:
                # fecha uma tira
                if len(strip) >= 3:
                    for k in range(len(strip) - 2):
                        i0, i1, i2 = strip[k], strip[k+1], strip[k+2]
                        # alterna winding
                        if k % 2 == 0:
                            a, b, c = i0, i1, i2
                        else:
                            a, b, c = i1, i0, i2

                        v0, v1, v2 = vtx(a), vtx(b), vtx(c)
                        p0 = GL._project_vertex(v0)
                        p1 = GL._project_vertex(v1)
                        p2 = GL._project_vertex(v2)
                        if p0 is None or p1 is None or p2 is None:
                            continue

                        e0, e1, e2 = to_eye(v0), to_eye(v1), to_eye(v2)
                        N = np.cross(e1 - e0, e2 - e0)
                        n = np.linalg.norm(N)
                        if n != 0.0:
                            N = (N / n).tolist()
                        else:
                            N = [0.0, 0.0, 1.0]

                        a0 = {"rgb": None, "uv": None, "n": N}
                        a1 = {"rgb": None, "uv": None, "n": N}
                        a2 = {"rgb": None, "uv": None, "n": N}
                        GL._raster_triangle_persp(p0, p1, p2, a0, a1, a2, colors, tex=None, tri_uv_area=None)
                strip = []
            else:
                strip.append(idx)


    @staticmethod
    def indexedFaceSet(coord, coordIndex, colorPerVertex, color, colorIndex,
                    texCoord, texCoordIndex, colors, current_texture):
        """
        Renderiza um IndexedFaceSet com sombreamento suave (normal por vértice)
        e correção de perspectiva. Suporta UV via texCoord/texCoordIndex.
        """

        # --- helpers -------------------------------------------------------------
        def to_eye(v, MV):
            P = np.array([v[0], v[1], v[2], 1.0], dtype=float)
            e = MV @ P
            return e[:3]

        def project(v):
            # Usa o projetor já existente do seu GL
            return GL._project_vertex(v)

        # --- prepara dados -------------------------------------------------------
        nverts = len(coord) // 3
        Vpos = [np.array(coord[i*3:i*3+3], dtype=float) for i in range(nverts)]

        # matriz mundo->câmera para normais
        V = getattr(GL, "_V", np.eye(4))
        M = getattr(GL, "_M", np.eye(4))
        MV = V @ M

        # lista de faces (índices) + início correspondente em texCoordIndex
        faces = []
        face = []
        poly_start_tc = 0
        for idx in coordIndex:
            if idx == -1:
                if len(face) >= 3:
                    faces.append((face[:], poly_start_tc))
                    if isinstance(texCoordIndex, (list, tuple)):
                        poly_start_tc += len(face) + 1  # avança no array de tcIndex
                face = []
            else:
                face.append(int(idx))
        # (se o arquivo não terminar com -1)
        if len(face) >= 3:
            faces.append((face[:], poly_start_tc))

        # --- acumula normal por vértice (no espaço da câmera) -------------------
        Vn = np.zeros((nverts, 3), dtype=float)
        for f, _tcstart in faces:
            # fan triangulation
            for k in range(1, len(f) - 1):
                i0, i1, i2 = f[0], f[k], f[k+1]
                e0 = to_eye(Vpos[i0], MV)
                e1 = to_eye(Vpos[i1], MV)
                e2 = to_eye(Vpos[i2], MV)
                n = np.cross(e1 - e0, e2 - e0)
                ln = np.linalg.norm(n)
                if ln > 0.0:
                    n /= ln
                # acumula nos 3 vértices
                Vn[i0] += n
                Vn[i1] += n
                Vn[i2] += n
        # normaliza
        for i in range(nverts):
            ln = np.linalg.norm(Vn[i])
            if ln > 0.0:
                Vn[i] /= ln
            else:
                Vn[i] = np.array([0.0, 0.0, 1.0], dtype=float)

        # --- util: pega UVs para um triângulo -----------------------------------
        def get_uv_for_corner(idx_vertex, pos_in_poly, tcstart):
            if texCoord is None:
                return None
            # Se há texCoordIndex, usa-o (por-canto)
            if isinstance(texCoordIndex, (list, tuple)) and len(texCoordIndex) > 0:
                ii = tcstart + pos_in_poly
                if 0 <= ii < len(texCoordIndex):
                    uvi = texCoordIndex[ii]
                    if uvi is not None and uvi >= 0 and len(texCoord) >= (uvi+1)*2:
                        return [float(texCoord[uvi*2]), float(texCoord[uvi*2+1])]
                return None
            # Caso contrário, assume que texCoord segue os índices de coord
            if len(texCoord) >= (idx_vertex+1)*2:
                return [float(texCoord[idx_vertex*2]), float(texCoord[idx_vertex*2+1])]
            return None

        # área em UV para mipmap (|u10 x u20|)
        def tri_uv_area(uv0, uv1, uv2):
            if (uv0 is None) or (uv1 is None) or (uv2 is None):
                return None
            u10 = np.array(uv1) - np.array(uv0)
            u20 = np.array(uv2) - np.array(uv0)
            return abs(u10[0]*u20[1] - u10[1]*u20[0])

        # textura (apenas se for dict com mipmaps)
        tex_obj = current_texture if (isinstance(current_texture, dict) and ("levels" in current_texture)) else None

        # --- rasteriza cada triângulo -------------------------------------------
        for f, tcstart in faces:
            for k in range(1, len(f) - 1):
                i0, i1, i2 = f[0], f[k], f[k+1]

                # projeção
                p0 = project(Vpos[i0]);  p1 = project(Vpos[i1]);  p2 = project(Vpos[i2])
                if (p0 is None) or (p1 is None) or (p2 is None):
                    continue

                # atributos por-vértice
                uv0 = get_uv_for_corner(i0, 0,  tcstart)
                uv1 = get_uv_for_corner(i1, k,  tcstart)
                uv2 = get_uv_for_corner(i2, k+1, tcstart)
                triA = tri_uv_area(uv0, uv1, uv2)

                a0 = {"rgb": None, "uv": uv0, "n": Vn[i0].tolist()}
                a1 = {"rgb": None, "uv": uv1, "n": Vn[i1].tolist()}
                a2 = {"rgb": None, "uv": uv2, "n": Vn[i2].tolist()}

                # chama o raster (com guarda de textura)
                GL._raster_triangle_persp(p0, p1, p2, a0, a1, a2, colors,
                                        tex=tex_obj, tri_uv_area=triA)




    @staticmethod
    def box(size, colors):
        """Função usada para renderizar Boxes."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Box

        print("Box : size = {0}".format(size)) # imprime no terminal pontos
        print("Box : colors = {0}".format(colors)) # imprime no terminal as cores

        gpu.GPU.draw_pixel([10, 10], gpu.GPU.RGB8, [255, 255, 255])  # altera pixel

    @staticmethod
    def sphere(radius, colors):
        """Função usada para renderizar Esferas."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Sphere

        print("Sphere : radius = {0}".format(radius)) # imprime no terminal o raio da esfera
        print("Sphere : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def cone(bottomRadius, height, colors):
        """Função usada para renderizar Cones."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cone

        print("Cone : bottomRadius = {0}".format(bottomRadius)) # imprime no terminal o raio da base do cone
        print("Cone : height = {0}".format(height)) # imprime no terminal a altura do cone
        print("Cone : colors = {0}".format(colors)) # imprime no terminal as cores

    @staticmethod
    def cylinder(radius, height, colors):
        """Função usada para renderizar Cilindros."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/geometry3D.html#Cylinder

        print("Cylinder : radius = {0}".format(radius)) # imprime no terminal o raio do cilindro
        print("Cylinder : height = {0}".format(height)) # imprime no terminal a altura do cilindro
        print("Cylinder : colors = {0}".format(colors)) # imprime no terminal as cores


    @staticmethod
    def navigationInfo(headlight):
        """Características físicas do avatar do visualizador e do modelo de visualização."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/navigation.html#NavigationInfo

        print("NavigationInfo : headlight = {0}".format(headlight)) # imprime no terminal
        GL._headlight_on = bool(headlight)  # (mantém só uma atribuição)

    @staticmethod
    def directionalLight(ambientIntensity, color, intensity, direction, **_kwargs):
        col = [1.0, 1.0, 1.0]
        if isinstance(color, (list, tuple)) and len(color) >= 3:
            col = [float(color[0]), float(color[1]), float(color[2])]

        I = 1.0 if intensity is None else float(intensity)
        amb = 0.0 if ambientIntensity is None else float(ambientIntensity)

        d = np.array(direction if direction is not None else [0.0, 0.0, -1.0], dtype=float)
        n = np.linalg.norm(d)
        d = d / n if n != 0.0 else np.array([0.0, 0.0, -1.0], dtype=float)

        GL._directional_light = {
            "ambientIntensity": amb,
            "intensity": I,
            "color": col,
            "direction_world": d,
        }

    @staticmethod
    def get_directional_light():
        return getattr(GL, "_directional_light", None)


    @staticmethod
    def pointLight(ambientIntensity, color, intensity, location):
        """Luz pontual."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/lighting.html#PointLight

        print("PointLight : ambientIntensity = {0}".format(ambientIntensity))
        print("PointLight : color = {0}".format(color)) # imprime no terminal
        print("PointLight : intensity = {0}".format(intensity)) # imprime no terminal
        print("PointLight : location = {0}".format(location)) # imprime no terminal

    @staticmethod
    def fog(visibilityRange, color):
        """Névoa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/environmentalEffects.html#Fog

        print("Fog : color = {0}".format(color)) # imprime no terminal
        print("Fog : visibilityRange = {0}".format(visibilityRange))

    @staticmethod
    def timeSensor(cycleInterval, loop):
        """Gera eventos conforme o tempo passa."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/time.html#TimeSensor

        print("TimeSensor : cycleInterval = {0}".format(cycleInterval)) # imprime no terminal
        print("TimeSensor : loop = {0}".format(loop))

        if cycleInterval is None or cycleInterval <= 0:
            return 0.0

        epoch = time.time()  # segundos desde a época
        if loop:
            fraction_changed = (epoch % cycleInterval) / cycleInterval
        else:
            if not hasattr(GL, "_ts_start"):
                GL._ts_start = epoch
            t = epoch - GL._ts_start
            fraction_changed = t / cycleInterval
            if fraction_changed >= 1.0:
                fraction_changed = 1.0
        return float(max(0.0, min(1.0, fraction_changed)))


    @staticmethod
    def splinePositionInterpolator(set_fraction, key, keyValue, closed):
        """Interpola não linearmente entre uma lista de vetores 3D."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#SplinePositionInterpolator

        print("SplinePositionInterpolator : set_fraction = {0}".format(set_fraction))
        print("SplinePositionInterpolator : key = {0}".format(key)) # imprime no terminal
        print("SplinePositionInterpolator : keyValue = {0}".format(keyValue))
        print("SplinePositionInterpolator : closed = {0}".format(closed))

        if not key or not keyValue or len(keyValue) < 3:
            return [0.0, 0.0, 0.0]

        K = list(map(float, key))
        P = [np.array(keyValue[i:i+3], dtype=float) for i in range(0, len(keyValue), 3)]
        n = len(K)

        if len(P) != n:
            t = float(max(0.0, min(1.0, set_fraction)))
            p0, p1 = P[0], P[-1]
            v = (1.0 - t) * p0 + t * p1
            return [float(v[0]), float(v[1]), float(v[2])]

        tglob = float(max(0.0, min(1.0, set_fraction)))

        if n == 1:
            v = P[0]
            return [float(v[0]), float(v[1]), float(v[2])]
        if n == 2:
            t = (tglob - K[0]) / max(1e-8, (K[1] - K[0]))
            t = max(0.0, min(1.0, t))
            v = (1.0 - t) * P[0] + t * P[1]
            return [float(v[0]), float(v[1]), float(v[2])]

        i = 0
        while i < n - 2 and tglob > K[i+1]:
            i += 1
        denom = max(1e-8, K[i+1] - K[i])
        t = (tglob - K[i]) / denom
        t = max(0.0, min(1.0, t))

        def idx(j):
            if closed:
                return (j + n) % n
            return max(0, min(n-1, j))

        p0 = P[idx(i-1)]
        p1 = P[idx(i  )]
        p2 = P[idx(i+1)]
        p3 = P[idx(i+2)]

        m1 = 0.5 * (p2 - p0)
        m2 = 0.5 * (p3 - p1)

        t2 = t * t
        t3 = t2 * t
        h00 =  2*t3 - 3*t2 + 1
        h10 =      t3 - 2*t2 + t
        h01 = -2*t3 + 3*t2
        h11 =      t3 -   t2

        v = h00*p1 + h10*m1 + h01*p2 + h11*m2
        return [float(v[0]), float(v[1]), float(v[2])]


    @staticmethod
    def orientationInterpolator(set_fraction, key, keyValue):
        """Interpola entre uma lista de valores de rotação especificos."""
        # https://www.web3d.org/specifications/X3Dv4/ISO-IEC19775-1v4-IS/Part01/components/interpolators.html#OrientationInterpolator

        print("OrientationInterpolator : set_fraction = {0}".format(set_fraction))
        print("OrientationInterpolator : key = {0}".format(key)) # imprime no terminal
        print("OrientationInterpolator : keyValue = {0}".format(keyValue))

        import math
        def aa_to_q(ax, ay, az, ang):
            v = np.array([ax, ay, az], dtype=float)
            n = np.linalg.norm(v)
            if n == 0.0:
                v = np.array([0.0, 0.0, 1.0], dtype=float)
                n = 1.0
            v = v / n
            h = 0.5 * float(ang)
            s = math.sin(h)
            return np.array([math.cos(h), v[0]*s, v[1]*s, v[2]*s], dtype=float)

        def q_normalize(q):
            return q / (np.linalg.norm(q) or 1.0)

        def q_to_aa(q):
            q = q_normalize(q)
            w, x, y, z = q
            w = max(-1.0, min(1.0, float(w)))
            ang = 2.0 * math.acos(w)
            s = math.sqrt(max(0.0, 1.0 - w*w))
            if s < 1e-8:
                axis = np.array([0.0, 0.0, 1.0], dtype=float)
            else:
                axis = np.array([x/s, y/s, z/s], dtype=float)
            return [float(axis[0]), float(axis[1]), float(axis[2]), float(ang)]

        def slerp(q0, q1, t):
            q0 = q_normalize(q0); q1 = q_normalize(q1)
            dot = float(np.dot(q0, q1))
            if dot < 0.0:
                q1 = -q1
                dot = -dot
            if dot > 0.9995:
                q = q0 + t*(q1 - q0)
                return q_normalize(q)
            theta0 = math.acos(max(-1.0, min(1.0, dot)))
            sin0 = math.sin(theta0)
            theta = theta0 * t
            s0 = math.sin(theta0 - theta) / sin0
            s1 = math.sin(theta) / sin0
            return s0*q0 + s1*q1

        if not key or not keyValue or len(keyValue) < 4:
            return [0, 0, 1, 0]

        K = list(map(float, key))
        R = [keyValue[i:i+4] for i in range(0, len(keyValue), 4)]
        n = len(K)
        if len(R) != n:
            r = R[0]
            return [float(r[0]), float(r[1]), float(r[2]), float(r[3])]

        tglob = float(max(0.0, min(1.0, set_fraction)))

        if n == 1 or tglob <= K[0]:
            r = R[0]
            return [float(r[0]), float(r[1]), float(r[2]), float(r[3])]
        if tglob >= K[-1]:
            r = R[-1]
            return [float(r[0]), float(r[1]), float(r[2]), float(r[3])]

        i = 0
        while i < n - 2 and tglob > K[i+1]:
            i += 1

        t = (tglob - K[i]) / max(1e-8, (K[i+1] - K[i]))
        t = max(0.0, min(1.0, t))

        q0 = aa_to_q(*R[i])
        q1 = aa_to_q(*R[i+1])

        q = slerp(q0, q1, t)
        ax, ay, az, ang = q_to_aa(q)
        value_changed = [ax, ay, az, ang]

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
        if not hasattr(GL, "_zbuf") or GL._zbuf.shape != (GL.height, GL.width):
            GL._zbuf = np.full((GL.height, GL.width), 1.0, dtype=float)

    @staticmethod
    def _clear_software_buffers():
        """Limpa z-buffer de software para 1.0 (chamar a cada frame)."""
        if hasattr(GL, "_zbuf"):
            GL._zbuf.fill(1.0)

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
        """
        Rasteriza triângulo com correção de perspectiva para z, cores, UV e NORMAL.
        Aceita material dict (diffuseColor/emissiveColor/transparency).
        'tex' pode ser None, ou um dict com {'levels': [...]} (mipmaps). Lista/str serão ignorados.
        """

        # --- Z-BUFFER ------------------------------------------------------------
        if not hasattr(GL, "_zbuf") or GL._zbuf.shape != (GL.height, GL.width):
            GL._zbuf = np.full((GL.height, GL.width), 1.0, dtype=float)

        # --- COR BASE DO MATERIAL ------------------------------------------------
        diff = material.get('diffuseColor', None) if isinstance(material, dict) else None
        emis = material.get('emissiveColor', None) if isinstance(material, dict) else None
        if diff is not None:
            base_rgb_mat = np.array([diff[0]*255.0, diff[1]*255.0, diff[2]*255.0], dtype=float)
        elif emis is not None:
            base_rgb_mat = np.array([emis[0]*255.0, emis[1]*255.0, emis[2]*255.0], dtype=float)
        else:
            base_rgb_mat = np.array([255.0, 255.0, 255.0], dtype=float)

        alpha_mat = 1.0 - float(material.get('transparency', 0.0)) if isinstance(material, dict) else 1.0

        # --- BBOX NA TELA --------------------------------------------------------
        minx = int(max(0, np.floor(min(p0["sx"], p1["sx"], p2["sx"]))))
        maxx = int(min(GL.width  - 1, np.ceil (max(p0["sx"], p1["sx"], p2["sx"]))))
        miny = int(max(0, np.floor(min(p0["sy"], p1["sy"], p2["sy"]))))
        maxy = int(min(GL.height - 1, np.ceil (max(p0["sy"], p1["sy"], p2["sy"]))))

        # área assinada (edge function)
        A = GL._edge((p0["sx"], p0["sy"]), (p1["sx"], p1["sy"]), (p2["sx"], p2["sy"]))
        if A == 0:
            return
        tri_screen_area = abs(A)

        # --- MIPMAP (apenas se 'tex' for dict com níveis) ------------------------
        level_img = None
        if isinstance(tex, dict) and ("levels" in tex) and (tri_uv_area is not None):
            L = GL._estimate_mip_level(tex, tri_uv_area, tri_screen_area)
            level_img = tex["levels"][L]

        # --- ATRIBUTOS POR-VÉRTICE ----------------------------------------------
        w0, w1, w2 = p0["invw"], p1["invw"], p2["invw"]
        z0, z1, z2 = p0["z_ndc"], p1["z_ndc"], p2["z_ndc"]

        # --- LOOP DE PIXELS ------------------------------------------------------
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

                # barycentr. não normalizada
                l0 = e0 / A; l1 = e1 / A; l2 = e2 / A

                # denom perspectiva
                denom = l0*w0 + l1*w1 + l2*w2
                if denom == 0:
                    continue

                # Z com correção de perspectiva
                z_ndc = (l0*z0*w0 + l1*z1*w1 + l2*z2*w2) / denom
                z01   = 0.5 * z_ndc + 0.5
                zcur  = GL._zbuf[y, x]

                # --- COR BASE (por-vértice opcional) -----------------------------
                if a0.get("rgb") is not None:
                    cr = (l0*a0["rgb"][0]*w0 + l1*a1["rgb"][0]*w1 + l2*a2["rgb"][0]*w2) / denom
                    cg = (l0*a0["rgb"][1]*w0 + l1*a1["rgb"][1]*w1 + l2*a2["rgb"][1]*w2) / denom
                    cb = (l0*a0["rgb"][2]*w0 + l1*a1["rgb"][2]*w1 + l2*a2["rgb"][2]*w2) / denom
                    rgb_base = np.array([cr, cg, cb], dtype=float)
                else:
                    rgb_base = base_rgb_mat.copy()

                rgb = rgb_base.copy()
                alpha_tex = 1.0

                # --- TEXTURA (se houver mip level válido + UVs) ------------------
                if (level_img is not None) and (a0.get("uv") is not None):
                    u = (l0*a0["uv"][0]*w0 + l1*a1["uv"][0]*w1 + l2*a2["uv"][0]*w2) / denom
                    v = (l0*a0["uv"][1]*w0 + l1*a1["uv"][1]*w1 + l2*a2["uv"][1]*w2) / denom
                    texel = GL._sample_tex_level(level_img, u, v)
                    if texel.shape[0] >= 3:
                        rgb = rgb_base * (texel[:3] / 255.0)
                    if texel.shape[0] == 4:
                        alpha_tex = float(texel[3]) / 255.0

                # --- NORMAL COM CORREÇÃO DE PERSPECTIVA --------------------------
                if a0.get("n") is not None:
                    Nx = (l0*a0["n"][0]*w0 + l1*a1["n"][0]*w1 + l2*a2["n"][0]*w2) / denom
                    Ny = (l0*a0["n"][1]*w0 + l1*a1["n"][1]*w1 + l2*a2["n"][1]*w2) / denom
                    Nz = (l0*a0["n"][2]*w0 + l1*a1["n"][2]*w1 + l2*a2["n"][2]*w2) / denom
                    Nn = np.array([Nx, Ny, Nz], dtype=float)
                    nlen = np.linalg.norm(Nn)
                    if nlen > 0.0:
                        Nn /= nlen

                    # luz direcional (em espaço da câmera)
                    light = getattr(GL, "_directional_light", None)
                    if light is not None:
                        # converte direção da luz mundo->olho
                        Vrot = getattr(GL, "_V", np.eye(4))[:3, :3]
                        Lw   = np.array(light["direction_world"], dtype=float)
                        Lv   = - (Vrot @ (Lw / (np.linalg.norm(Lw) or 1.0)))  # para a superfície
                        ndotl = max(0.0, float(np.dot(Nn, Lv)))

                        ambI = float(material.get('ambientIntensity', 0.2)) if isinstance(material, dict) else 0.2
                        ambL = float(light.get('ambientIntensity', 0.0))
                        I    = float(light.get('intensity', 1.0))
                        Lcol = np.array(light['color'], dtype=float)  # 0..1

                        total = np.clip(ambI + ambL + I*ndotl, 0.0, 1.0) * Lcol
                        rgb = rgb * total  # rgb está em 0..255

                # --- ALPHA -------------------------------------------------------
                alpha = float(np.clip(alpha_mat * alpha_tex, 0.0, 1.0))

                # --- DEPTH/BLEND -------------------------------------------------
                if alpha >= 0.999:
                    if z01 < zcur:
                        gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8,
                                        [int(np.clip(rgb[0],0,255)),
                                            int(np.clip(rgb[1],0,255)),
                                            int(np.clip(rgb[2],0,255))])
                        GL._zbuf[y, x] = z01
                else:
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
        dx = abs(x1 - x0)
        dy = abs(y1 - y0)
        
        if dx > dy:
            if x0 > x1:
                x0, x1 = x1, x0
                y0, y1 = y1, y0
            
            for x in range(x0, x1 + 1):
                y = int(y0 + (y1 - y0) * (x - x0) / (x1 - x0))
                if 0 <= x < GL.width and 0 <= y < GL.height:
                    gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, color)
        else:
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
        vertices = [(x1, y1), (x2, y2), (x3, y3)]
        vertices.sort(key=lambda v: v[1])
        
        x1, y1 = vertices[0]
        x2, y2 = vertices[1]
        x3, y3 = vertices[2]
        
        for y in range(int(y1), int(y3) + 1):
            if 0 <= y < GL.height:
                if y <= y2:
                    if y2 != y1:
                        x_start = int(x1 + (x2 - x1) * (y - y1) / (y2 - y1))
                        x_end = int(x1 + (x3 - x1) * (y - y1) / (y3 - y1))
                    else:
                        x_start = x1
                        x_end = x1
                else:
                    if y3 != y2:
                        x_start = int(x2 + (x3 - x2) * (y - y2) / (y3 - y2))
                        x_end = int(x1 + (x3 - x1) * (y - y1) / (y3 - y1))
                    else:
                        x_start = x2
                        x_end = x2
                
                if x_start > x_end:
                    x_start, x_end = x_end, x_start
                
                for x in range(x_start, x_end + 1):
                    if 0 <= x < GL.width:
                        gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, color)



    # [ADICIONADO] ---------- Iluminação (difusa + especular) flat por triângulo ----------
    @staticmethod
    def _apply_lighting_flat(world_center, world_normal, material):
        # material
        diff = np.array(material.get('diffuseColor', [1.0, 1.0, 1.0]), dtype=float)
        spec = np.array(material.get('specularColor', [0.0, 0.0, 0.0]), dtype=float)
        ambI = float(material.get('ambientIntensity', 0.0))
        shin = float(material.get('shininess', 0.2))  # X3D: 0..1
        shin_exp = max(1.0, 1.0 + shin * 127.0)

        N = np.array(world_normal, dtype=float)
        N = N / (np.linalg.norm(N) or 1.0)

        cam_pos = getattr(GL, "_cam_pos", np.array([0.0, 0.0, 10.0], dtype=float))
        V = cam_pos - np.array(world_center, dtype=float)
        V = V / (np.linalg.norm(V) or 1.0)

        Li = np.zeros(3, dtype=float)

        if getattr(GL, "_headlight_on", False):
            dH = getattr(GL, "_cam_fwd", np.array([0.0, 0.0, -1.0], dtype=float))
            Lvec = -dH / (np.linalg.norm(dH) or 1.0)
            ndotl = max(0.0, float(np.dot(N, Lvec)))
            R = 2.0 * ndotl * N - Lvec
            rdotv = max(0.0, float(np.dot(R, V))) ** shin_exp
            Li += (ndotl * diff) + (rdotv * spec)

        if hasattr(GL, "_lights"):
            for L in GL._lights:
                if L.get("type") != "dir":
                    continue
                col = L["color"] * float(L["intensity"])
                Ldir = -L["dir_world"]
                ndotl = max(0.0, float(np.dot(N, Ldir)))
                R = 2.0 * ndotl * N - Ldir
                rdotv = max(0.0, float(np.dot(R, V))) ** shin_exp
                ambient = L["ambientIntensity"] * diff
                Li += ambient * col + (ndotl * diff + rdotv * spec) * col

        return np.clip(Li, 0.0, 1.0)

    @staticmethod
    def begin_frame(bg=(0, 0, 0)):
        """Zera o framebuffer e o z-buffer no início de cada frame."""
        # --- limpa cor (tenta API da sua GPU; se não existir, faz fallback lento) ---
        ok = False
        try:
            # muitas versões têm isso:
            gpu.GPU.clear(list(map(int, bg)))
            ok = True
        except Exception:
            try:
                # algumas versões usam clear_color(r,g,b)
                gpu.GPU.clear_color(int(bg[0]), int(bg[1]), int(bg[2]))
                ok = True
            except Exception:
                pass

        if not ok:
            # fallback: pinta pixel a pixel (funciona em qualquer caso)
            rr, gg, bb = int(bg[0]), int(bg[1]), int(bg[2])
            for y in range(GL.height):
                for x in range(GL.width):
                    gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, [rr, gg, bb])

        # --- z-buffer = 1.0 (longe) ---
        if not hasattr(GL, "_zbuf") or GL._zbuf.shape != (GL.height, GL.width):
            GL._zbuf = np.full((GL.height, GL.width), 1.0, dtype=float)
        else:
            GL._zbuf.fill(1.0)


    def render(self):
        GL.begin_frame((0, 0, 0))  # fundo preto
        # ... depois configura câmera/luz e desenha a cena normalmente ...
        self.scene.render()
        return gpu.GPU.get_image()  # ou o que a sua interface espera retornar

