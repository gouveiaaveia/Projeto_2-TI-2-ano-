# Author: Marco Simoes
# Adapted from Java's implementation of Rui Pedro Paiva
# Teoria da Informacao, LEI, 2022

import sys
import numpy as np
from huffmantree import HuffmanTree


class GZIPHeader:

    ID1 = ID2 = CM = FLG = XFL = OS = 0
    MTIME = []
    lenMTIME = 4
    mTime = 0

    # bits 0, 1, 2, 3 and 4, respectively (remaining 3 bits: reserved)
    FLG_FTEXT = FLG_FHCRC = FLG_FEXTRA = FLG_FNAME = FLG_FCOMMENT = 0

    # FLG_FTEXT --> ignored (usually 0)
    # if FLG_FEXTRA == 1
    XLEN, extraField = [], []
    lenXLEN = 2

    # if FLG_FNAME == 1
    fName = ''  # ends when a byte with value 0 is read

    # if FLG_FCOMMENT == 1
    fComment = ''  # ends when a byte with value 0 is read

    # if FLG_HCRC == 1
    HCRC = []

    def read(self, f):
        ''' reads and processes the Huffman header from file. Returns 0 if no error, -1 otherwise '''

        # ID 1 and 2: fixed values
        self.ID1 = f.read(1)[0]
        if self.ID1 != 0x1f: return -1  # error in the header

        self.ID2 = f.read(1)[0]
        if self.ID2 != 0x8b: return -1  # error in the header

        # CM - Compression Method: must be the value 8 for deflate
        self.CM = f.read(1)[0]
        if self.CM != 0x08: return -1  # error in the header

        # Flags
        self.FLG = f.read(1)[0]

        # MTIME
        self.MTIME = [0] * self.lenMTIME
        self.mTime = 0
        for i in range(self.lenMTIME):
            self.MTIME[i] = f.read(1)[0]
            self.mTime += self.MTIME[i] << (8 * i)

        # XFL (not processed...)
        self.XFL = f.read(1)[0]

        # OS (not processed...)
        self.OS = f.read(1)[0]

        # --- Check Flags
        self.FLG_FTEXT = self.FLG & 0x01
        self.FLG_FHCRC = (self.FLG & 0x02) >> 1
        self.FLG_FEXTRA = (self.FLG & 0x04) >> 2
        self.FLG_FNAME = (self.FLG & 0x08) >> 3
        self.FLG_FCOMMENT = (self.FLG & 0x10) >> 4

        # FLG_EXTRA
        if self.FLG_FEXTRA == 1:
            # read 2 bytes XLEN + XLEN bytes de extra field
            # 1st byte: LSB, 2nd: MSB
            self.XLEN = [0] * self.lenXLEN
            self.XLEN[0] = f.read(1)[0]
            self.XLEN[1] = f.read(1)[0]
            self.xlen = self.XLEN[1] << 8 + self.XLEN[0]

            # read extraField and ignore its values
            self.extraField = f.read(self.xlen)

        def read_str_until_0(f):
            s = ''
            while True:
                c = f.read(1)[0]
                if c == 0:
                    return s
                s += chr(c)

        # FLG_FNAME
        if self.FLG_FNAME == 1:
            self.fName = read_str_until_0(f)

        # FLG_FCOMMENT
        if self.FLG_FCOMMENT == 1:
            self.fComment = read_str_until_0(f)

        # FLG_FHCRC (not processed...)
        if self.FLG_FHCRC == 1:
            self.HCRC = f.read(2)

        return 0


class GZIP:
    ''' class for GZIP decompressing file (if compressed with deflate) '''

    gzh = None
    gzFile = ''
    fileSize = origFileSize = -1
    numBlocks = 0
    f = None

    bits_buffer = 0
    available_bits = 0

    def __init__(self, filename):
        self.gzFile = filename
        self.f = open(filename, 'rb')
        self.f.seek(0, 2)
        self.fileSize = self.f.tell()
        self.f.seek(0)

    def decompress(self):
        ''' main function for decompressing the gzip file with deflate algorithm '''

        numBlocks = 0

        # get original file size: size of file before compression
        origFileSize = self.getOrigFileSize()
        print(origFileSize)

        # read GZIP header
        error = self.getHeader()
        if error != 0:
            print('Formato invalido!')
            return

        # show filename read from GZIP header
        print(self.gzh.fName)

        # MAIN LOOP - decode block by block
        BFINAL = 0
        while not BFINAL == 1:

            BFINAL = self.readBits(1)

            BTYPE = self.readBits(2)
            if BTYPE != 2:
                print('Error: Block %d not coded with Huffman Dynamic coding' % (numBlocks + 1))
                return

            # --- STUDENTS --- ADD CODE HERE

            # read HLIT, HDIST and HCLEN
            HLIT = self.readBits(5) + 257
            HDIST = self.readBits(5) + 1
            HCLEN = self.readBits(4) + 4

            # criar uma lista com 19 elementos, todos a 0
            lista_valores_comprimento_HCLEN = np.zeros(19, dtype=int)
            ordem_da_lista = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]

            # ler os valores de HCLEN
            for i in range(HCLEN):
                lista_valores_comprimento_HCLEN[ordem_da_lista[i]] = self.readBits(3)

            # Listas de valores decimais e binários (strings)
            lista_valores_decimal_comprimentos = np.full(19, '', dtype='U30')
            lista_valores_binario_comprimentos = np.full(19, '', dtype='U30')

            # Lista única de comprimentos e o valor inicial
            lista_comprimentos_comprimentos = np.unique(lista_valores_comprimento_HCLEN[lista_valores_comprimento_HCLEN > 0])

            print(f"\nLista de valores de comprimento: {lista_valores_comprimento_HCLEN}")
            #print(f"Lista de comprimentos: {lista_comprimentos}")

            # Atribuir valores decimais e binários
            self.valores_decimais(lista_comprimentos_comprimentos, lista_valores_decimal_comprimentos, lista_valores_comprimento_HCLEN)

            # Converter para binários com base no comprimento de cada elemento
            self.valores_binarios(lista_valores_decimal_comprimentos, lista_valores_binario_comprimentos, lista_valores_comprimento_HCLEN)
            print(f"\nLista de valores binários: {lista_valores_binario_comprimentos}")

            #criar a árvore de Huffman (semana 2)
            hft = HuffmanTree()
            verbose = False
            self.arvore_huffman(lista_valores_binario_comprimentos, hft, verbose)

            #semana 3
            #ponto 4
            lista_comprimentos_HLIT = np.zeros(286, dtype=int)
            
            self.leitura_comprimentos_codigo(lista_comprimentos_HLIT, hft, HLIT)
            print(f"\nLista de comprimentos HLIT: {lista_comprimentos_HLIT}")

            lista_valores_decimal_lit = np.full(286, '', dtype='U30')
            lista_valores_binario_lit = np.full(286, '', dtype='U30')

            lista_comprimentos_lit = np.unique(lista_comprimentos_HLIT[lista_comprimentos_HLIT > 0])

            self.valores_decimais(lista_comprimentos_lit, lista_valores_decimal_lit, lista_comprimentos_HLIT)

            self.valores_binarios(lista_valores_decimal_lit, lista_valores_binario_lit, lista_comprimentos_HLIT)
            print(f"\nLista de valores binários lit: {lista_valores_binario_lit}")

            hft2 = HuffmanTree() # arvore dos literais
            verbose = False
            self.arvore_huffman(lista_valores_binario_lit, hft2, verbose)

            #ponto 5
            lista_comprimentos_HDIST = np.zeros(30, dtype=int)
            self.leitura_comprimentos_codigo(lista_comprimentos_HDIST, hft, HDIST)
            print(f"\nLista de comprimentos HDIST: {lista_comprimentos_HDIST}")

            lista_valores_decimal_dist = np.full(30, '', dtype='U30')
            lista_valores_binario_dist = np.full(30, '', dtype='U30')

            lista_comprimentos_dist = np.unique(lista_comprimentos_HDIST[lista_comprimentos_HDIST > 0])

            self.valores_decimais(lista_comprimentos_dist, lista_valores_decimal_dist, lista_comprimentos_HDIST)
        
            self.valores_binarios(lista_valores_decimal_dist, lista_valores_binario_dist, lista_comprimentos_HDIST)
            print(f"\nLista de valores binários dist: {lista_valores_binario_dist}")

            hft3 = HuffmanTree() # arvore das distâncias
            verbose = False
            self.arvore_huffman(lista_valores_binario_dist, hft3, verbose)

            #ponto 7
            max_dist = 32768  # Distância máxima permitida

            # Processa a descompactação
            lista_descompactada = self.descompactacao(hft2, hft3, max_dist)

            self.escrever_em_ficheiro(lista_descompactada, fileName[0:-3])

            # update number of blocks read
            numBlocks += 1

        # close file

        self.f.close()
        print("End: %d block(s) analyzed." % numBlocks)


    def valores_decimais(self, lista_comprimentos, lista_valores_decimal, lista_valores_comprimento):
        valor = 0
        comprimento_anterior = lista_comprimentos[0]

        # Atribuir valores decimais e binários
        for comprimento in lista_comprimentos:
            # Verificar diferença entre comprimentos
            diferenca = comprimento - comprimento_anterior

            if diferenca > 1:
                # Incrementar e multiplicar o valor por 2 para cada diferença maior que 1
                for _ in range(diferenca - 1):
                    valor = (valor) << 1

            # Posições onde o comprimento atual está presente
            indices = np.where(lista_valores_comprimento == comprimento)[0]
            
            # Atribuir valores decimais sequenciais para os índices
            lista_valores_decimal[indices] = [str(valor + i) for i in range(len(indices))]

            # Atualizar o valor base e o comprimento anterior
            valor += len(indices)
            valor <<= 1  # Ajustar para o próximo comprimento
            comprimento_anterior = comprimento


    def valores_binarios(self, lista_valores_decimal, lista_valores_binario, lista_valores_comprimento):
        '''Converte valores decimais para binários com o preenchimento adequado de zeros à esquerda'''
        for i, val in enumerate(lista_valores_decimal):
            if val:
                # Converte o valor decimal para binário
                binario = bin(int(val))[2:]  # Remove o prefixo '0b'
            
                # Preenche com zeros à esquerda para que o comprimento binário seja igual ao comprimento especificado
                print(f"{lista_valores_comprimento[i]} {binario}")
                lista_valores_binario[i] = binario.zfill(lista_valores_comprimento[i]) 



    def arvore_huffman(self, lista_valores_binario, hft, verbose):
        for code in range(len(lista_valores_binario)):
                if lista_valores_binario[code] != '':
                    hft.addNode(lista_valores_binario[code], code, verbose)


    def leitura_comprimentos_codigo(self, lista_comprimento, hft, tamanho):
            hft.resetCurNode()
            i = 0
            posicao_array = 0
            code = ""
            while posicao_array < tamanho:
                leitura = self.readBits(1)
                code += str(leitura)
                pos = hft.nextNode(str(leitura))

                if pos == -1:
                    print("Code '" + code + "' not found!!!")
                    break
                elif pos == -2:
                    continue
                else:
                    if pos == 16:
                        valor = self.readBits(2) + 3
                        lista_comprimento[posicao_array:posicao_array + valor] = lista_comprimento[posicao_array - 1]
                        posicao_array += valor
                    elif pos == 17:
                        valor = self.readBits(3) + 3
                        lista_comprimento[posicao_array:posicao_array + valor] = 0
                        posicao_array += valor
                    elif pos == 18:
                        valor = self.readBits(7) + 11
                        lista_comprimento[posicao_array:posicao_array + valor] = 0
                        posicao_array += valor
                    else:
                        lista_comprimento[posicao_array] = pos
                        posicao_array += 1
                    code = ""
                    hft.resetCurNode()


    def descompactacao(self, hft_literais, hft_distancias, max_dist):

        comprimentos = {257: (0, 3), 258: (0, 4), 259: (0, 5), 260: (0, 6), 261: (0, 7), 262: (0, 8), 263: (0, 9),
                        264: (0, 10), 265: (1, 11), 266: (1, 13), 267: (1, 15), 268: (1, 17), 269: (2, 19), 270: (2, 23), 
                        271: (2, 27), 272: (2, 31), 273: (3, 35), 274: (3, 43), 275: (3, 51), 276: (3, 59), 
                        277: (4, 67), 278: (4, 83), 279: (4, 99), 280: (4, 115), 281: (5, 131), 282: (5, 163), 
                        283: (5, 195), 284: (5, 227), 285: (0, 258)}

        distanicas = {0: (0, 1), 1: (0, 2), 2: (0, 3), 3: (0, 4), 4: (1, 5), 5: (1, 7), 6: (2, 9), 7: (2, 13), 
                    8: (3, 17), 9: (3, 25), 10: (4, 33), 11: (4, 49), 12: (5, 65), 13: (5, 97), 14: (6, 129), 
                    15: (6, 193), 16: (7, 257), 17: (7, 385), 18: (8, 513), 19: (8, 769), 20: (9, 1025), 
                    21: (9, 1537), 22: (10, 2049), 23: (10, 3073), 24: (11, 4097), 25: (11, 6145), 26: (12, 8193), 
                    27: (12, 12289), 28: (13, 16385), 29: (13, 24577)}

        lista_descompactada = []
        hft_literais.resetCurNode()
        code = ""
        comprimento = 0

        while True:
            leitura = self.readBits(1)
            code += str(leitura)
            pos = hft_literais.nextNode(str(leitura))

            if pos == -1:
                print("Code " + code + " not found!!")
                break
            elif pos == -2:
                continue

            hft_literais.resetCurNode()

            if pos == 256:
                print("End of block")
                break

            if pos < 256:
                lista_descompactada.append(pos)
                hft_literais.resetCurNode()
                while len(lista_descompactada) > max_dist:
                    del lista_descompactada[0]
                continue

            else:
                comprimento = self.readBits(comprimentos[pos][0]) if comprimentos[pos][0] > 0 else 0
                comprimento += comprimentos[pos][1]

                hft_distancias.resetCurNode()

                while True:
                    leitura2 = self.readBits(1)
                    pos2 = hft_distancias.nextNode(str(leitura2))

                    if pos2 == -1:
                        print("Code " + str(leitura2) + " not found!!!")
                        break
                    elif pos2 == -2:
                        continue
                    else:
                        dist = self.readBits(distanicas[pos2][0]) if distanicas[pos2][0] > 0 else 0 
                        dist += distanicas[pos2][1]
                        
                        for i in range(comprimento):
                            if len(lista_descompactada) - dist >= 0:
                                valor = lista_descompactada[len(lista_descompactada) - dist]
                                lista_descompactada.append(valor)
                                while len(lista_descompactada) > max_dist:
                                    del lista_descompactada[0]
                        break

        return lista_descompactada

    def escrever_em_ficheiro(self, lista_descompactada, ficheiro_saida):
        with open(ficheiro_saida, "ab+") as file:
            file.write(bytearray(lista_descompactada))


    def getOrigFileSize(self):
        ''' reads file size of original file (before compression) - ISIZE '''

        # saves current position of file pointer
        fp = self.f.tell()

        # jumps to end-4 position
        self.f.seek(self.fileSize - 4)

        # reads the last 4 bytes (LITTLE ENDIAN)
        sz = 0
        for i in range(4):
            sz += self.f.read(1)[0] << (8 * i)

        # restores file pointer to its original position
        self.f.seek(fp)

        return sz

    def getHeader(self):
        ''' reads GZIP header'''

        self.gzh = GZIPHeader()
        header_error = self.gzh.read(self.f)
        return header_error

    def readBits(self, n, keep=False):
        ''' reads n bits from bits_buffer. if keep = True, leaves bits in the buffer for future accesses '''

        while n > self.available_bits:
            self.bits_buffer = self.f.read(1)[0] << self.available_bits | self.bits_buffer
            self.available_bits += 8

        mask = (2 ** n) - 1
        value = self.bits_buffer & mask

        if not keep:
            self.bits_buffer >>= n
            self.available_bits -= n

        return value


if __name__ == '__main__':

    # gets filename from command line if provided
    fileName = "sample_audio.mp3.gz"
    if len(sys.argv) > 1:
        fileName = sys.argv[1]

    # decompress file
    gz = GZIP(fileName)
    gz.decompress()
