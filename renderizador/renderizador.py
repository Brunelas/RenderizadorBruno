#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
Renderizador X3D.

Desenvolvido por: Luciano Soares <lpsoares@insper.edu.br>
Disciplina: Computação Gráfica
Data: 28 de Agosto de 2020
"""

import os           # Para rotinas do sistema operacional
import argparse     # Para tratar os parâmetros da linha de comando

import gl           # Recupera rotinas de suporte ao X3D

import interface    # Janela de visualização baseada no Matplotlib
import gpu          # Simula os recursos de uma GPU

import x3d          # Faz a leitura do arquivo X3D, gera o grafo de cena e faz traversal
import scenegraph   # Imprime o grafo de cena no console

LARGURA = 60  # Valor padrão para largura da tela
ALTURA = 40   # Valor padrão para altura da tela


class Renderizador:
    """Realiza a renderização da cena informada."""

    def __init__(self):
        """Definindo valores padrão."""
        self.width = LARGURA
        self.height = ALTURA
        self.x3d_file = ""
        self.image_file = "tela.png"
        self.scene = None
        self.framebuffers = {}

    def setup(self):
        """Configura o sistema para a renderização."""
        # Configurando color buffers para exibição na tela

        # ADICIONADO: fator de supersampling (2x2)
        self.ssaa = 2  # pode deixar fixo em 2 (enunciado pede 2x2 por padrão)

        # Cria uma (2) posição de FrameBuffer na GPU
        fbo = gpu.GPU.gen_framebuffers(2) # ADICIONADO: pedir 2 FBOs (um para tela, outro para SSAA)

        # Define o atributo FRONT como o FrameBuffe principal
        self.framebuffers["FRONT"] = fbo[0]

        # ADICIONADO: framebuffer “SSAA” onde realmente vamos rasterizar em resolução maior
        self.framebuffers["SSAA"] = fbo[1]

        # Define que a posição criada será usada para desenho e leitura
        gpu.GPU.bind_framebuffer(gpu.GPU.FRAMEBUFFER, self.framebuffers["FRONT"])
        # Opções:
        # - DRAW_FRAMEBUFFER: Faz o bind só para escrever no framebuffer
        # - READ_FRAMEBUFFER: Faz o bind só para leitura no framebuffer
        # - FRAMEBUFFER: Faz o bind para leitura e escrita no framebuffer

        # Aloca memória no FrameBuffer para um tipo e tamanho especificado de buffer

        # Memória de Framebuffer para canal de cores
        gpu.GPU.framebuffer_storage(
            self.framebuffers["FRONT"],
            gpu.GPU.COLOR_ATTACHMENT,
            gpu.GPU.RGB8,
            self.width,
            self.height
        )

        # Descomente as seguintes linhas se for usar um Framebuffer para profundidade
        gpu.GPU.framebuffer_storage(
            self.framebuffers["FRONT"],
            gpu.GPU.DEPTH_ATTACHMENT,
            gpu.GPU.DEPTH_COMPONENT32F,
            self.width,
            self.height
        )

        # ADICIONADO: aloca o framebuffer “SSAA” com 2x a resolução na cor e na profundidade
        ss_w, ss_h = self.width * self.ssaa, self.height * self.ssaa
        gpu.GPU.bind_framebuffer(gpu.GPU.FRAMEBUFFER, self.framebuffers["SSAA"])
        gpu.GPU.framebuffer_storage(
            self.framebuffers["SSAA"],
            gpu.GPU.COLOR_ATTACHMENT,
            gpu.GPU.RGB8,
            ss_w,
            ss_h
        )
        gpu.GPU.framebuffer_storage(
            self.framebuffers["SSAA"],
            gpu.GPU.DEPTH_ATTACHMENT,
            gpu.GPU.DEPTH_COMPONENT32F,
            ss_w,
            ss_h
        )
    
        # Opções:
        # - COLOR_ATTACHMENT: alocações para as cores da imagem renderizada
        # - DEPTH_ATTACHMENT: alocações para as profundidades da imagem renderizada
        # Obs: Você pode chamar duas vezes a rotina com cada tipo de buffer.

        # Tipos de dados:
        # - RGB8: Para canais de cores (Vermelho, Verde, Azul) 8bits cada (0-255)
        # - RGBA8: Para canais de cores (Vermelho, Verde, Azul, Transparência) 8bits cada (0-255)
        # - DEPTH_COMPONENT16: Para canal de Profundidade de 16bits (half-precision) (0-65535)
        # - DEPTH_COMPONENT32F: Para canal de Profundidade de 32bits (single-precision) (float)

        # Define cor que ira apagar o FrameBuffer quando clear_buffer() invocado
        gpu.GPU.clear_color([0, 0, 0])

        # Define a profundidade que ira apagar o FrameBuffer quando clear_buffer() invocado
        # Assuma 1.0 o mais afastado e -1.0 o mais próximo da camera
        gpu.GPU.clear_depth(1.0)

        # ADICIONADO: vamos desenhar no SSAA (alta resolução)
        #             bind para leitura+escrita no FBO SSAA
        gpu.GPU.bind_framebuffer(gpu.GPU.FRAMEBUFFER, self.framebuffers["SSAA"])

        # Definindo tamanho do Viewport para renderização
        #self.scene.viewport(width=self.width, height=self.height)

        # ADICIONADO: ajustar viewport e GL para a resolução SSAA
        self.scene.viewport(width=ss_w, height=ss_h)
        # o GL usa width/height para mapear NDC->tela; atualizamos para ssaa
        import gl
        gl.GL.width = ss_w
        gl.GL.height = ss_h

        gl.GL._alloc_software_buffers()  # [ADICIONADO] aloca z-buffer de software na resolução SSAA
        



    def pre(self):
        """Rotinas pré renderização."""
        # Função invocada antes do processo de renderização iniciar.

        # Limpa o frame buffers atual
        gpu.GPU.clear_buffer()

        gl.GL._clear_software_buffers()  # [ADICIONADO] z-buffer de software volta para +inf

        if hasattr(gl.GL, "_zbuf"):
            gl.GL._zbuf.fill(1.0)  # 1.0 = mais longe (compatível com NDC)
        # Recursos que podem ser úteis:
        # Define o valor do pixel no framebuffer: draw_pixel(coord, mode, data)
        # Retorna o valor do pixel no framebuffer: read_pixel(coord, mode)

    def pos(self):
        """Rotinas pós renderização."""
        # Função invocada após o processo de renderização terminar.

        # Essa é uma chamada conveniente para manipulação de buffers
        # ao final da renderização de um frame. Como por exemplo, executar
        # downscaling da imagem.

        # ADICIONADO: downscale 2x2 do SSAA -> FRONT
        if hasattr(self, "ssaa") and self.ssaa == 2 and "SSAA" in self.framebuffers:
            w, h = self.width, self.height
            ss_w, ss_h = w * 2, h * 2

            # Bind leitura no SSAA e escrita no FRONT
            gpu.GPU.bind_framebuffer(gpu.GPU.READ_FRAMEBUFFER, self.framebuffers["SSAA"])
            gpu.GPU.bind_framebuffer(gpu.GPU.DRAW_FRAMEBUFFER, self.framebuffers["FRONT"])

            # varre cada pixel de saída (w x h) e faz média de um bloco 2x2 do SSAA
            for y in range(h):
                sy = 2 * y
                for x in range(w):
                    sx = 2 * x
                    # lê 4 amostras do SSAA
                    c00 = gpu.GPU.read_pixel([sx,     sy    ], gpu.GPU.RGB8)
                    c10 = gpu.GPU.read_pixel([sx + 1, sy    ], gpu.GPU.RGB8)
                    c01 = gpu.GPU.read_pixel([sx,     sy + 1], gpu.GPU.RGB8)
                    c11 = gpu.GPU.read_pixel([sx + 1, sy + 1], gpu.GPU.RGB8)
                    # média simples (box filter 2x2)
                    r = (c00[0] + c10[0] + c01[0] + c11[0]) // 4
                    g = (c00[1] + c10[1] + c01[1] + c11[1]) // 4
                    b = (c00[2] + c10[2] + c01[2] + c11[2]) // 4
                    gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, [int(r), int(g), int(b)])

        # (opcional) se quiser limpar o SSAA depois, pode fazê-lo aqui.

        # Método para a troca dos buffers (NÃO IMPLEMENTADO)
        # Esse método será utilizado na fase de implementação de animações
        gpu.GPU.swap_buffers()

    def mapping(self):
        """Mapeamento de funções para as rotinas de renderização."""
        # Rotinas encapsuladas na classe GL (Graphics Library)
        x3d.X3D.renderer["Polypoint2D"] = gl.GL.polypoint2D
        x3d.X3D.renderer["Polyline2D"] = gl.GL.polyline2D
        x3d.X3D.renderer["Circle2D"] = gl.GL.circle2D
        x3d.X3D.renderer["TriangleSet2D"] = gl.GL.triangleSet2D
        x3d.X3D.renderer["TriangleSet"] = gl.GL.triangleSet
        x3d.X3D.renderer["Viewpoint"] = gl.GL.viewpoint
        x3d.X3D.renderer["Transform_in"] = gl.GL.transform_in
        x3d.X3D.renderer["Transform_out"] = gl.GL.transform_out
        x3d.X3D.renderer["TriangleStripSet"] = gl.GL.triangleStripSet
        x3d.X3D.renderer["IndexedTriangleStripSet"] = gl.GL.indexedTriangleStripSet
        x3d.X3D.renderer["IndexedFaceSet"] = gl.GL.indexedFaceSet
        x3d.X3D.renderer["Box"] = gl.GL.box
        x3d.X3D.renderer["Sphere"] = gl.GL.sphere
        x3d.X3D.renderer["Cone"] = gl.GL.cone
        x3d.X3D.renderer["Cylinder"] = gl.GL.cylinder
        x3d.X3D.renderer["NavigationInfo"] = gl.GL.navigationInfo
        x3d.X3D.renderer["DirectionalLight"] = gl.GL.directionalLight
        x3d.X3D.renderer["PointLight"] = gl.GL.pointLight
        x3d.X3D.renderer["Fog"] = gl.GL.fog
        x3d.X3D.renderer["TimeSensor"] = gl.GL.timeSensor
        x3d.X3D.renderer["SplinePositionInterpolator"] = gl.GL.splinePositionInterpolator
        x3d.X3D.renderer["OrientationInterpolator"] = gl.GL.orientationInterpolator

    def render(self):
        """Laço principal de renderização."""
        self.pre()  # executa rotina pré renderização
        self.scene.render()  # faz o traversal no grafo de cena
        self.pos()  # executa rotina pós renderização
        return gpu.GPU.get_frame_buffer()

    def main(self):
        """Executa a renderização."""
        # Tratando entrada de parâmetro
        parser = argparse.ArgumentParser(add_help=False)   # parser para linha de comando
        parser.add_argument("-i", "--input", help="arquivo X3D de entrada")
        parser.add_argument("-o", "--output", help="arquivo 2D de saída (imagem)")
        parser.add_argument("-w", "--width", help="resolução horizonta", type=int)
        parser.add_argument("-h", "--height", help="resolução vertical", type=int)
        parser.add_argument("-g", "--graph", help="imprime o grafo de cena", action='store_true')
        parser.add_argument("-p", "--pause", help="começa simulação em pausa", action='store_true')
        parser.add_argument("-q", "--quiet", help="não exibe janela", action='store_true')
        args = parser.parse_args() # parse the arguments
        if args.input:
            self.x3d_file = args.input
        if args.output:
            self.image_file = args.output
        if args.width:
            self.width = args.width
        if args.height:
            self.height = args.height

        path = os.path.dirname(os.path.abspath(self.x3d_file))

        # Iniciando simulação de GPU
        gpu.GPU(self.image_file, path)

        # Abre arquivo X3D
        self.scene = x3d.X3D(self.x3d_file)

        # Iniciando Biblioteca Gráfica
        gl.GL.setup(
            self.width,
            self.height,
            near=0.01,
            far=1000
        )

        # Funções que irão fazer o rendering
        self.mapping()

        # Se no modo silencioso não configurar janela de visualização
        if not args.quiet:
            window = interface.Interface(self.width, self.height, self.x3d_file)
            self.scene.set_preview(window)

        # carrega os dados do grafo de cena
        if self.scene:
            self.scene.parse()
            if args.graph:
                scenegraph.Graph(self.scene.root)

        # Configura o sistema para a renderização.
        self.setup()

        # Se no modo silencioso salvar imagem e não mostrar janela de visualização
        if args.quiet:
            gpu.GPU.save_image()  # Salva imagem em arquivo
        else:
            window.set_saver(gpu.GPU.save_image)  # pasa a função para salvar imagens
            window.preview(args.pause, self.render)  # mostra visualização

if __name__ == '__main__':
    renderizador = Renderizador()
    renderizador.main()
