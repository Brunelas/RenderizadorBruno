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
import numpy as np


import gl           # Recupera rotinas de suporte ao X3D
from gl import GL

import interface    # Janela de visualização baseada no Matplotlib
import gpu          # Simula os recursos de uma GPU

import x3d          # Faz a leitura do arquivo X3D, gera o grafo de cena e faz traversal
import scenegraph   # Imprime o grafo de cena no console

LARGURA = 60  # Valor padrão para largura da tela
ALTURA = 40   # Valor padrão para altura da tela


class Renderizador:
    """Realiza a renderização da cena informada."""

    def __init__(self):
        self.width = LARGURA
        self.height = ALTURA
        self.x3d_file = ""
        self.image_file = "tela.png"
        self.scene = None
        self.framebuffers = {}
        self.ssaa = 1        # valor inicial; o setup() vai forçar 2
        # NÃO calcule ssaa_width/ssaa_height aqui; vamos calcular no setup()

    def setup(self):
        """Configura o sistema para a renderização."""

        # 2x SSAA (padrão do enunciado)
        self.ssaa = 2
        self.ssaa_width  = self.width  * self.ssaa
        self.ssaa_height = self.height * self.ssaa

        # Cria 2 FBOs: FRONT (final) e SSAA (alta resolução)
        fbo = gpu.GPU.gen_framebuffers(2)
        self.framebuffers["FRONT"] = fbo[0]
        self.framebuffers["SSAA"]  = fbo[1]

        # FRONT (tamanho da janela)
        gpu.GPU.bind_framebuffer(gpu.GPU.FRAMEBUFFER, self.framebuffers["FRONT"])
        gpu.GPU.framebuffer_storage(self.framebuffers["FRONT"], gpu.GPU.COLOR_ATTACHMENT, gpu.GPU.RGB8, self.width, self.height)
        gpu.GPU.framebuffer_storage(self.framebuffers["FRONT"], gpu.GPU.DEPTH_ATTACHMENT, gpu.GPU.DEPTH_COMPONENT32F, self.width, self.height)

        # SSAA (tamanho 2× da janela)
        gpu.GPU.bind_framebuffer(gpu.GPU.FRAMEBUFFER, self.framebuffers["SSAA"])
        gpu.GPU.framebuffer_storage(self.framebuffers["SSAA"], gpu.GPU.COLOR_ATTACHMENT, gpu.GPU.RGB8, self.ssaa_width, self.ssaa_height)
        gpu.GPU.framebuffer_storage(self.framebuffers["SSAA"], gpu.GPU.DEPTH_ATTACHMENT, gpu.GPU.DEPTH_COMPONENT32F, self.ssaa_width, self.ssaa_height)

        # clear padrão
        gpu.GPU.clear_color([0, 0, 0])
        gpu.GPU.clear_depth(1.0)

        # Vamos rasterizar no FBO SSAA
        gpu.GPU.bind_framebuffer(gpu.GPU.FRAMEBUFFER, self.framebuffers["SSAA"])

        # Ajusta GL e viewport para a resolução SSAA (ESSENCIAL!)
        gl.GL.setup(self.ssaa_width, self.ssaa_height, near=gl.GL.near, far=gl.GL.far)
        self.scene.viewport(width=self.ssaa_width, height=self.ssaa_height)

        # z-buffer de software no mesmo tamanho do raster SSAA
        gl.GL._alloc_software_buffers()
        



    def pre(self):
        """Rotinas pré renderização."""
        # garantir que estamos desenhando no FBO SSAA
        gpu.GPU.bind_framebuffer(gpu.GPU.FRAMEBUFFER, self.framebuffers["SSAA"])

        # limpa color+depth do FBO atual
        gpu.GPU.clear_buffer()

        # limpa z-buffer de software (0..1, onde 1.0 = mais longe)
        gl.GL._alloc_software_buffers()
        gl.GL._clear_software_buffers()

    def pos(self):
        """Rotinas pós renderização."""
        if self.ssaa == 2 and "SSAA" in self.framebuffers:
            w, h = self.width, self.height

            # lê do SSAA, escreve no FRONT
            gpu.GPU.bind_framebuffer(gpu.GPU.READ_FRAMEBUFFER, self.framebuffers["SSAA"])
            gpu.GPU.bind_framebuffer(gpu.GPU.DRAW_FRAMEBUFFER, self.framebuffers["FRONT"])

            for y in range(h):
                sy = 2 * y
                for x in range(w):
                    sx = 2 * x
                    c00 = gpu.GPU.read_pixel([sx,     sy    ], gpu.GPU.RGB8)
                    c10 = gpu.GPU.read_pixel([sx + 1, sy    ], gpu.GPU.RGB8)
                    c01 = gpu.GPU.read_pixel([sx,     sy + 1], gpu.GPU.RGB8)
                    c11 = gpu.GPU.read_pixel([sx + 1, sy + 1], gpu.GPU.RGB8)
                    # dentro de pos(), onde faz a média 2x2
                    r = int(round((c00[0] + c10[0] + c01[0] + c11[0]) / 4))
                    g = int(round((c00[1] + c10[1] + c01[1] + c11[1]) / 4))
                    b = int(round((c00[2] + c10[2] + c01[2] + c11[2]) / 4))
                    gpu.GPU.draw_pixel([x, y], gpu.GPU.RGB8, [int(r), int(g), int(b)])

        # deixe o FRONT como framebuffer “corrente”
        gpu.GPU.bind_framebuffer(gpu.GPU.FRAMEBUFFER, self.framebuffers["FRONT"])
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
        GL.begin_frame((0, 0, 0))  
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
        # gl.GL.setup(
        #     self.ssaa_width, self.ssaa_height, near=gl.GL.near, far=gl.GL.far
        # )

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
