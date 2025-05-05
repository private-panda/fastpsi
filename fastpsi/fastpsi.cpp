
#include "libOTe/Base/BaseOT.h"
#include "libOTe/TwoChooseOne/Kos/KosOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Kos/KosOtExtSender.h"
#include "libOTe/TwoChooseOne/KosDot/KosDotExtReceiver.h"
#include "libOTe/TwoChooseOne/KosDot/KosDotExtSender.h"
#include "libOTe/TwoChooseOne/Iknp/IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Iknp/IknpOtExtSender.h"


#include "libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Silent/SilentOtExtSender.h"

#include "util.h"
#include "coproto/Socket/AsioSocket.h"

#include "../volePSI/Paxos.h"

using namespace osuCrypto;
using namespace volePSI;


void runReceiver(int nx, int ny, int omega, cp::AsioSocket &chl, bool verbose=0)
{

	oc::block setSeed1(toBlock(54321)), setSeed2(sysRandomSeed());
	std::vector<block> receiverSet(nx);
	PRNG prng1_set(setSeed1), prng2_set(setSeed2);
	int intersectionSize = 100;
	prng1_set.get<block>(receiverSet.data(), intersectionSize);
	prng2_set.get<block>(receiverSet.data()+intersectionSize, nx-intersectionSize);

	int columnBlockNumber = (omega+127)/128;
	// std::vector<block> senderValues(n);

	Timer timer;
	auto s = timer.setTimePoint("start");

	//encode to to get okvs table;
	int w = 3;
	int ssp = 40;
	auto dt = PaxosParam::Binary;

	Baxos paxos;
	paxos.init(nx, 1<<15, w, ssp, dt, block(0, 0));

	std::cout<<"The expansion rate: "<<(paxos.size()*1.0/nx)<<std::endl;
	oc::Matrix<block> valC(nx, columnBlockNumber), pax(paxos.size(), columnBlockNumber);
	// std::vector<std::vector<block>> valC(n, std::vector<block>(4)), pax(pp.size(), std::vector<block>(4));
	int okvsTableSize = paxos.size();
	// oc::PRNG codewordPrng(toBlock(11111));
	// for(int i=0; i<nx; i++){
	// 	codewordPrng.SetSeed(receiverSet[i]);
	// 	codewordPrng.get<block>(valC[i]);
	// }
	std::vector<block> valC_tmp(nx);
	std::vector<oc::AES> codewordsGenerator(columnBlockNumber);
	oc::block baseSeed(toBlock(11111));

	for(int i=0; i<codewordsGenerator.size(); i++){
		codewordsGenerator[i].setKey(baseSeed+toBlock(i));
		codewordsGenerator[i].ecbEncBlocks(receiverSet, valC_tmp);
		for(int j=0; j<nx; j++){
			valC[j][i] = valC_tmp[j];
		}
	}

	paxos.solve<block>(receiverSet, valC, pax, nullptr, dt);

	timer.setTimePoint("offline");

	PRNG prng(toBlock(12345));

	SilentOtExtSender  sender;

	DefaultBaseOT base;
	BitVector bv(128);
	std::array<block, 128> baseMsg;
	bv.randomize(prng);

	// perform the base To, call sync_wait to block until they have completed.
	cp::sync_wait(base.receive(bv, baseMsg, prng, chl));
	sender.setBaseOts(baseMsg, bv);

	//should be divided by 128 for transpose
	u64 numCols = (omega+127)/128*128;

	// construct a vector to stored the random send messages.
	std::vector<std::array<block, 2>> sMsgs(numCols);
	// perform the OTs and write the random OTs to msgs.
	cp::sync_wait(sender.send(sMsgs, prng, chl));


	auto e = timer.setTimePoint("ot finish");
	auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(e - s).count();

	auto com = (chl.bytesSent() + chl.bytesReceived())*1.0/(1<<20);

	lout << " n=" << Color::Green << numCols << " " << milli << " ms  " << com << " MB" << std::endl << Color::Default;
	std::cout<<"The superBlkSize: "<<superBlkSize<<std::endl;

	int rowBlockNumber = (paxos.size()+127)/128;
	//expand T0=sMsgs[0] and T1=sMsgs[1] to 

	static const u64 superBlkSize(8);
	int numOtExt = paxos.size();

	// this will be used as temporary buffers of 128 columns,
	// each containing 1024 bits. Once transposed, they will be copied
	// into the T1, T0 buffers for long term storage.
	std::array<std::array<block, superBlkSize>, 128> t0;
	std::array<std::array<block, superBlkSize>, 128> t1;

	// we are going to process OTs in blocks of 128 * superblkSize messages.
	u64 numSuperBlocks = ((numOtExt + 127) / 128 + superBlkSize - 1) / superBlkSize;

	// We need two matrices, T0 and T1. These will hold the expanded and transposed
	// rows that we got the using the base OTs as PRNG seed.
	oc::Matrix<block> mT0(numOtExt, numCols / 128);
	oc::Matrix<block> mT1(numOtExt, numCols / 128);

	// The is the index of the last correction value u = T0 ^ T1 ^ c(w)
	// that was sent to the sender.
	int mCorrectionIdx = 0;

	std::vector<std::array<AES, 2>> mGens(numCols);
	for(int i=0; i<mGens.size(); i++){
		mGens[i][0].setKey(sMsgs[i][0]);
		mGens[i][1].setKey(sMsgs[i][1]);
	}
	std::vector<u64> mGensBlkIdx(numCols, 0);

	// the index of the OT that has been completed.
	u64 doneIdx = 0;

	// NOTE: We do not transpose a bit-matrix of size numCol * numCol.
	//   Instead we break it down into smaller chunks. We do 128 columns
	//   times 8 * 128 rows at a time, where 8 = superBlkSize. This is done for
	//   performance reasons. The reason for 8 is that most CPUs have 8 AES vector
	//   lanes, and so its more efficient to encrypt (aka prng) 8 blocks at a time.
	//   So that's what we do.
	for (u64 superBlkIdx = 0; superBlkIdx < numSuperBlocks; ++superBlkIdx)
	{
		// compute at what row does the user want us to stop.
		// The code will still compute the transpose for these
		// extra rows, but it is thrown away.
		u64 stopIdx
			= doneIdx
			+ std::min<u64>(u64(128) * superBlkSize, numOtExt - doneIdx);


		for (u64 i = 0; i < numCols / 128; ++i)
		{

			for (u64 tIdx = 0, colIdx = i * 128; tIdx < 128; ++tIdx, ++colIdx)
			{
				// generate the column indexed by colIdx. This is done with
				// AES in counter mode acting as a PRNG. We don't use the normal
				// PRNG interface because that would result in a data copy when
				// we move it into the T0,T1 matrices. Instead we do it directly.
				mGens[colIdx][0].ecbEncCounterMode(mGensBlkIdx[colIdx], superBlkSize, ((block*)t0.data() + superBlkSize * tIdx));
				mGens[colIdx][1].ecbEncCounterMode(mGensBlkIdx[colIdx], superBlkSize, ((block*)t1.data() + superBlkSize * tIdx));

				// increment the counter mode idx.
				mGensBlkIdx[colIdx] += superBlkSize;
			}

			// transpose our 128 columns of 1024 bits. We will have 1024 rows,
			// each 128 bits wide.
			transpose128x1024(t0);
			transpose128x1024(t1);

			// This is the index of where we will store the matrix long term.
			// doneIdx is the starting row. i is the offset into the blocks of 128 bits.
			// __restrict isn't crucial, it just tells the compiler that this pointer
			// is unique and it shouldn't worry about pointer aliasing.
			block* __restrict mT0Iter = mT0.data() + mT0.stride() * doneIdx + i;
			block* __restrict mT1Iter = mT1.data() + mT1.stride() * doneIdx + i;

			for (u64 rowIdx = doneIdx, j = 0; rowIdx < stopIdx; ++j)
			{
				// because we transposed 1024 rows, the indexing gets a bit weird. But this
				// is the location of the next row that we want. Keep in mind that we had long
				// **contiguous** columns.
				block* __restrict t0Iter = ((block*)t0.data()) + j;
				block* __restrict t1Iter = ((block*)t1.data()) + j;

				// do the copy!
				for (u64 k = 0; rowIdx < stopIdx && k < 128; ++rowIdx, ++k)
				{
					*mT0Iter = *(t0Iter);
					*mT1Iter = *(t1Iter);

					t0Iter += superBlkSize;
					t1Iter += superBlkSize;

					mT0Iter += mT0.stride();
					mT1Iter += mT0.stride();
				}
			}
		}

		doneIdx = stopIdx;
	}
	timer.setTimePoint("transpose T0 and T1");
	// std::cout<<"finishing transposing T0 and T1"<<std::endl;

	int columnBytesNumber = (omega+7)/8;
	//send the correction to the sender
	std::vector<u8> sendBuffer(numOtExt*columnBytesNumber);
	u8* sendBuffer_ptr = sendBuffer.data();
	auto mT1_ptr = mT1.data();
	for(int i=0; i<numOtExt; i++){
		for(int j=0; j<columnBlockNumber; j++){
			mT1[i][j] ^= (mT0[i][j] ^ pax[i][j]);
		}
		std::memcpy(sendBuffer_ptr, (u8 *)(mT1_ptr), columnBytesNumber);
		sendBuffer_ptr += columnBytesNumber;
		mT1_ptr += mT1.stride();
	}
	cp::sync_wait(chl.send(sendBuffer));


	// std::cout<<"Receiver transpose the paxos table "<<std::endl;
	timer.setTimePoint("send the correction");
	double comSent_oprf = (chl.bytesSent())*1.0/(1<<20);
	double comRecv_oprf = (chl.bytesReceived())*1.0/(1<<20);
	lout << "Communication (MB): sent & recvd & total " << Color::Green << comSent_oprf << " " << comRecv_oprf << " " << (comSent_oprf+comRecv_oprf) << std::endl << Color::Default;

	// lout << "Communication (MB): sent & recvd " << Color::Green << (chl.bytesSent())*1.0/(1<<20) << " " << (chl.bytesReceived())*1.0/(1<<20) << std::endl << Color::Default;

	//send the correction t0 ^ sMsgs[0] ^ sMsgs[1]
	oc::Matrix<block> vals(nx, columnBlockNumber);
	// paxos.template decode<block>(receiverSet, vals, mT0);
	paxos.decode<block>(receiverSet, vals, mT0, dt);

	// std::cout<<"Receiver transpose the paxos Q"<<std::endl;
	timer.setTimePoint("get the T0 values");

	// std::cout<<"The vals size: "<<vals.size()<<", "<<vals[0].size()<<std::endl;
	//computet the oprf values and put it into a map
	std::unordered_map<u32, std::pair<block,int>> oprfSet;
	// std::unordered_map<block, int> oprfSet;
	int outputByteNumber = ((int)(ssp + log2(nx) + log2(ny))+7)/8;
	RandomOracle H(outputByteNumber);
	// std::vector<block> tmpResult(columnBlockNumber);
	block result;

	auto vals_ptr = vals.data();
	for(int i=0; i<nx; i++){
		H.Reset();
		H.Update((u8 *)(vals_ptr), columnBytesNumber);
		vals_ptr += vals.stride();
		H.Final((u8 *)&result);

		oprfSet[*(u32 *)&result] = std::pair<oc::block, int>(result, i);
		// oprfSet[result] = i;
	}

	// std::cout<<"Receiver gets the oprf values"<<std::endl;
	timer.setTimePoint("oprf");

	// std::cout<<"The outputByteNumber: "<<outputByteNumber<<std::endl;
	std::vector<u8> recvBuffer(ny*outputByteNumber);
	// std::vector<block> recvBuffer((n*outputByteNumber+15)/16);
	// std::cout<<"The buffer size: "<<recvBuffer.size()<<std::endl;
	cp::sync_wait(chl.recv(recvBuffer));

	u8 *recvPtr = recvBuffer.data();
	std::vector<block> intersection;
	// oc::block oprf_value;
	for(int i=0; i<ny; i++){
		auto iter = oprfSet.find(*(uint32_t *)recvPtr);
		if(iter != oprfSet.end() && std::memcmp(recvPtr, (u8 *)&(iter->second.first), outputByteNumber) == 0){
			intersection.emplace_back(receiverSet[iter->second.second]);
		}
		
		recvPtr += outputByteNumber;
	}
	std::cout<<"The intersection size: "<<intersection.size()<<std::endl;

	timer.setTimePoint("intersection");

	double comSent = (chl.bytesSent())*1.0/(1<<20);
	double comRecv = (chl.bytesReceived())*1.0/(1<<20);

	std::cout<<"Receiving sender's OPRF value (MB): "<<(comRecv - comRecv_oprf)<<std::endl;

	com = (chl.bytesSent() + chl.bytesReceived())*1.0/(1<<20);
	lout << " n=" << Color::Green << omega << " " << milli << " ms  " << com << " MB" << std::endl << Color::Default;
	lout << "Communication (MB): sent & recvd & total " << Color::Green << comSent << " " << comRecv << " " << (comSent+comRecv)<<std::endl << Color::Default;

	lout << " **** receiver ****\n" << timer << std::endl;


	chl.close();
}

void runSender(int nx, int ny, int omega, cp::AsioSocket &chl, bool verbose=0)
{
	oc::block setSeed1(toBlock(54321)), setSeed2(sysRandomSeed());
	std::vector<block> senderSet(ny);
	PRNG prng1_set(setSeed1), prng2_set(setSeed2);
	int intersectionSize = 100;
	prng1_set.get<block>(senderSet.data(), intersectionSize);
	prng2_set.get<block>(senderSet.data()+intersectionSize, ny-intersectionSize);

	Timer timer;
	auto s = timer.setTimePoint("start");

	int columnBlockNumber = (omega+127)/128;

	//encode to to get okvs table;
	int w = 3;
	int ssp = 40;
	auto dt = PaxosParam::Binary;

	Baxos paxos;
	paxos.init(nx, 1<<15, w, ssp, dt, block(0, 0));

	auto cols = 0;

	// PaxosParam pp(nx, w, ssp, dt);
	std::cout<<"The expansion rate: "<<(paxos.size()*1.0/nx)<<std::endl;
	oc::Matrix<block> valC(ny, columnBlockNumber), pax(paxos.size(), columnBlockNumber);

	std::vector<block> valC_tmp(ny);
	std::vector<oc::AES> codewordsGenerator(columnBlockNumber);
	oc::block baseSeed(toBlock(11111));

	for(int i=0; i<codewordsGenerator.size(); i++){
		codewordsGenerator[i].setKey(baseSeed+toBlock(i));
		codewordsGenerator[i].ecbEncBlocks(senderSet, valC_tmp);
		for(int j=0; j<ny; j++){
			valC[j][i] = valC_tmp[j];
		}
	}

	// PRNG prng(sysRandomSeed());
	PRNG prng(toBlock(12345));

	SilentOtExtReceiver  receiver;

	// Now compute the base OTs, we need to set them on the first pair of extenders.
	// In real code you would only have a sender or reciever, not both. But we do
	// here just showing the example.

	DefaultBaseOT base;
	std::array<std::array<block, 2>, 128> baseMsg;

	// perform the base To, call sync_wait to block until they have completed.
	cp::sync_wait(base.send(baseMsg, prng, chl));
	receiver.setBaseOts(baseMsg);

	u64 numCols = (omega+127)/128*128;

	// construct the choices that we want.
	BitVector choice(numCols);
	// in this case pick random messages.
	choice.randomize(prng);

	// construct a vector to stored the received messages.
	std::vector<block> rMsgs(numCols);
	// perform  totalOTs random OTs, the results will be written to msgs.
	cp::sync_wait(receiver.receive(choice, rMsgs, prng, chl));

	auto e = timer.setTimePoint("finish");
	auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(e - s).count();

	auto com = (chl.bytesSent() + chl.bytesReceived())*1.0/(1<<20);

	// lout << " number of columns = " << Color::Green << numCols << " " << milli << " ms  " << com << " MB" << std::endl << Color::Default;

	int okvsTableSize = paxos.size();
	int rowBlockNumber = (okvsTableSize+127)/128;
	// std::vector<block> expand_msgs(rowBlockNumber);

	//the number of column should be divided by 128

	// round up
	int numOTExt = okvsTableSize;
	// numOTExt = ((numOTExt + 127) / 128) * 128;
	std::cout<<"The table size: "<<numOTExt<<std::endl;

	std::vector<AES> mGens(numCols);
	for(int i=0; i<numCols; i++){
		mGens[i].setKey(rMsgs[i]);
	}

	// We need two matrices, one for the senders matrix T^i_{b_i} and
	// one to hold the the correction values. This is sometimes called
	// the u = T0 + T1 + C matrix in the papers.
	oc::Matrix<block> mT(numOTExt, numCols/ 128);
	//char c;
	//chl.recv(&c, 1);

	oc::Matrix<block> mCorrectionVals(numOTExt, numCols / 128);

	// The receiver will send us correction values, this is the index of
	// the next one they will send.
	int mCorrectionIdx = 0;

	// we are going to process OTs in blocks of 128 * superblkSize messages.
	u64 numSuperBlocks = ((numOTExt + 127) / 128 + superBlkSize - 1) / superBlkSize;

	std::vector<u64> mGensBlkIdx(numCols, 0);

	// the index of the last OT that we have completed.
	u64 doneIdx = 0;

	// a temp that will be used to transpose the sender's matrix
	std::array<std::array<block, superBlkSize>, 128> t;

	for (u64 superBlkIdx = 0; superBlkIdx < numSuperBlocks; ++superBlkIdx)
	{
		// compute at what row does the user want use to stop.
		// the code will still compute the transpose for these
		// extra rows, but it is thrown away.
		u64 stopIdx
			= doneIdx
			+ std::min<u64>(u64(128) * superBlkSize, mT.bounds()[0] - doneIdx);

		// transpose 128 columns at at time. Each column will be 128 * superBlkSize = 1024 bits long.
		for (u64 i = 0; i < numCols / 128; ++i)
		{
			// generate the columns using AES-NI in counter mode.
			for (u64 tIdx = 0, colIdx = i * 128; tIdx < 128; ++tIdx, ++colIdx)
			{
				mGens[colIdx].ecbEncCounterMode(mGensBlkIdx[colIdx], superBlkSize, ((block*)t.data() + superBlkSize * tIdx));
				mGensBlkIdx[colIdx] += superBlkSize;
			}

			// transpose our 128 columns of 1024 bits. We will have 1024 rows,
			// each 128 bits wide.
			transpose128x1024(t);

			// This is the index of where we will store the matrix long term.
			// doneIdx is the starting row. i is the offset into the blocks of 128 bits.
			// __restrict isn't crucial, it just tells the compiler that this pointer
			// is unique and it shouldn't worry about pointer aliasing.
			block* __restrict mTIter = mT.data() + doneIdx * mT.stride() + i;

			for (u64 rowIdx = doneIdx, j = 0; rowIdx < stopIdx; ++j)
			{
				// because we transposed 1024 rows, the indexing gets a bit weird. But this
				// is the location of the next row that we want. Keep in mind that we had long
				// **contiguous** columns.
				block* __restrict tIter = (((block*)t.data()) + j);

				// do the copy!
				for (u64 k = 0; rowIdx < stopIdx && k < 128; ++rowIdx, ++k)
				{
					*mTIter = *tIter;

					tIter += superBlkSize;
					mTIter += mT.stride();
				}
			}

		}

		doneIdx = stopIdx;
	}

	timer.setTimePoint("transpose Q");
	// std::cout<<"finishing transposing Q"<<std::endl;

	int columnByteNumber = (omega+7)/8;
	std::vector<u8> recvBuffer(columnByteNumber*numOTExt);
	cp::sync_wait(chl.recv(recvBuffer));

	std::vector<block> choiceBlocks(columnBlockNumber);
	std::memcpy((u8 *)choiceBlocks.data(), (u8 *)choice.data(), columnByteNumber);

	//now receive the corrections
	auto mCorrection_ptr = mCorrectionVals.data();
	auto recvBuffter_ptr = recvBuffer.data();
	for(int i=0; i<numOTExt; i++){
		// cp::sync_wait(chl.recv(recvBuffer));
		// std::memcpy((u8 *)(mCorrectionVals.data()+mCorrectionVals.stride()*i), recvBuffer.data()+columnByteNumber*i, columnByteNumber);
		std::memcpy((u8 *)(mCorrection_ptr), recvBuffter_ptr, columnByteNumber);
		mCorrection_ptr += mCorrectionVals.stride();
		recvBuffter_ptr += columnByteNumber;

		for(int j=0; j<columnBlockNumber; j++){
			mT[i][j] ^= (mCorrectionVals[i][j] & choiceBlocks[j]);
		}
	}

	// std::cout<<"Sender gets the correction"<<std::endl;
	timer.setTimePoint("get correction");

	oc::Matrix<block> vals1(ny, columnBlockNumber);
	paxos.decode<block>(senderSet, vals1, mT, dt);

	int outputByteNumber = ((int)(ssp + log2(nx) + log2(ny))+7)/8;
	// std::cout<<"outputByteNumber: "<<outputByteNumber<<std::endl;
	RandomOracle H(outputByteNumber);
	std::vector<block> tmpResult(columnBlockNumber);
	std::vector<u8> oprfSendBuffer(outputByteNumber*ny);
	// std::vector<block> oprfSendBuffer((n*outputByteNumber+15)/16);

	// std::cout<<"The buffer size: "<<oprfSendBuffer.size()<<std::endl;

	std::vector<int> shuffle_indices(ny);
	for(int i=0; i<ny; i++){
		shuffle_indices[i] = i;
	}
	std::random_shuffle(shuffle_indices.begin(), shuffle_indices.end());

	u8 *oprtPtr = oprfSendBuffer.data();
	for(int i=0; i<ny; i++){
		H.Reset();
		for(int j=0; j<columnBlockNumber; j++){
			// tmpResult[j] = vals1[i][j] ^ ((vals2[i][j]^valC[i][j]) & choiceBlocks[j]);
			tmpResult[j] = vals1[i][j] ^ (valC[i][j] & choiceBlocks[j]);
			// H.Update(tmpResult[j]);
		}

		H.Update((u8 *)(tmpResult.data()), columnByteNumber);
		H.Final(oprtPtr+shuffle_indices[i]*outputByteNumber);

	}
	cp::sync_wait(chl.send(oprfSendBuffer));

	timer.setTimePoint("get oprf values");

	// double comSent = (chl.bytesSent())*1.0/(1<<20);
	// double comRecv = (chl.bytesReceived())*1.0/(1<<20);

	// com = (chl.bytesSent() + chl.bytesReceived())*1.0/(1<<20);
	// lout << " n=" << Color::Green << omega << " " << milli << " ms  " << com << " MB" << std::endl << Color::Default;
	// lout << "Communication (MB): sent & recvd " << Color::Green << comSent << " " << comRecv << std::endl << Color::Default;

	if (verbose)
	{
		lout << " **** sender ****\n" << timer << std::endl;
	};
	chl.close();
}


void runReceiverOpprf(int nx, int ny, int valueBlocks, int omega, cp::AsioSocket &chl, bool verbose=0)
{

	oc::block setSeed1(toBlock(54321)), setSeed2(sysRandomSeed());
	std::vector<block> receiverSet(nx);
	PRNG prng1_set(setSeed1), prng2_set(setSeed2);
	int intersectionSize = 100;
	prng1_set.get<block>(receiverSet.data(), intersectionSize);
	prng2_set.get<block>(receiverSet.data()+intersectionSize, nx-intersectionSize);

	int columnBlockNumber = (omega+127)/128;
	// std::vector<block> senderValues(n);

	Timer timer;
	auto s = timer.setTimePoint("start");

	//encode to to get okvs table;
	int w = 3;
	int ssp = 40;
	auto dt = PaxosParam::Binary;

	Baxos paxos;
	paxos.init(nx, 1<<15, w, ssp, dt, block(0, 0));

	std::cout<<"The expansion rate: "<<(paxos.size()*1.0/nx)<<std::endl;
	oc::Matrix<block> valC(nx, columnBlockNumber), pax(paxos.size(), columnBlockNumber);
	int okvsTableSize = paxos.size();

	std::vector<block> valC_tmp(nx);
	std::vector<oc::AES> codewordsGenerator(columnBlockNumber);
	oc::block baseSeed(toBlock(11111));

	for(int i=0; i<codewordsGenerator.size(); i++){
		codewordsGenerator[i].setKey(baseSeed+toBlock(i));
		codewordsGenerator[i].ecbEncBlocks(receiverSet, valC_tmp);
		for(int j=0; j<nx; j++){
			valC[j][i] = valC_tmp[j];
		}
	}

	paxos.solve<block>(receiverSet, valC, pax, nullptr, dt);

	timer.setTimePoint("offline");

	// PRNG prng(sysRandomSeed());
	PRNG prng(toBlock(12345));

	SilentOtExtSender  sender;
	// Now compute the base OTs, we need to set them on the first pair of extenders.
	// In real code you would only have a sender or reciever, not both. But we do
	// here just showing the example.
	DefaultBaseOT base;
	BitVector bv(128);
	std::array<block, 128> baseMsg;
	bv.randomize(prng);

	// perform the base To, call sync_wait to block until they have completed.
	cp::sync_wait(base.receive(bv, baseMsg, prng, chl));
	sender.setBaseOts(baseMsg, bv);

	//should be divided by 128 for transpose
	u64 numCols = (omega+127)/128*128;

	// construct a vector to stored the random send messages.
	std::vector<std::array<block, 2>> sMsgs(numCols);
	// perform the OTs and write the random OTs to msgs.
	cp::sync_wait(sender.send(sMsgs, prng, chl));


	auto e = timer.setTimePoint("ot finish");
	auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(e - s).count();

	auto com = (chl.bytesSent() + chl.bytesReceived())*1.0/(1<<20);

	lout << " nx=" <<nx<<", omega"<< Color::Green << omega << " " << milli << " ms  " << com << " MB" << std::endl << Color::Default;
	std::cout<<"The superBlkSize: "<<superBlkSize<<std::endl;

	int rowBlockNumber = (paxos.size()+127)/128;
	//expand T0=sMsgs[0] and T1=sMsgs[1] to 

	static const u64 superBlkSize(8);
	int numOtExt = paxos.size();

	// this will be used as temporary buffers of 128 columns,
	// each containing 1024 bits. Once transposed, they will be copied
	// into the T1, T0 buffers for long term storage.
	std::array<std::array<block, superBlkSize>, 128> t0;
	std::array<std::array<block, superBlkSize>, 128> t1;

	// we are going to process OTs in blocks of 128 * superblkSize messages.
	u64 numSuperBlocks = ((numOtExt + 127) / 128 + superBlkSize - 1) / superBlkSize;

	// We need two matrices, T0 and T1. These will hold the expanded and transposed
	// rows that we got the using the base OTs as PRNG seed.
	oc::Matrix<block> mT0(numOtExt, numCols / 128);
	oc::Matrix<block> mT1(numOtExt, numCols / 128);

	// The is the index of the last correction value u = T0 ^ T1 ^ c(w)
	// that was sent to the sender.
	int mCorrectionIdx = 0;

	std::vector<std::array<AES, 2>> mGens(numCols);
	for(int i=0; i<mGens.size(); i++){
		mGens[i][0].setKey(sMsgs[i][0]);
		mGens[i][1].setKey(sMsgs[i][1]);
	}
	std::vector<u64> mGensBlkIdx(numCols, 0);

	// the index of the OT that has been completed.
	u64 doneIdx = 0;

	// NOTE: We do not transpose a bit-matrix of size numCol * numCol.
	//   Instead we break it down into smaller chunks. We do 128 columns
	//   times 8 * 128 rows at a time, where 8 = superBlkSize. This is done for
	//   performance reasons. The reason for 8 is that most CPUs have 8 AES vector
	//   lanes, and so its more efficient to encrypt (aka prng) 8 blocks at a time.
	//   So that's what we do.
	for (u64 superBlkIdx = 0; superBlkIdx < numSuperBlocks; ++superBlkIdx)
	{
		// compute at what row does the user want us to stop.
		// The code will still compute the transpose for these
		// extra rows, but it is thrown away.
		u64 stopIdx
			= doneIdx
			+ std::min<u64>(u64(128) * superBlkSize, numOtExt - doneIdx);


		for (u64 i = 0; i < numCols / 128; ++i)
		{

			for (u64 tIdx = 0, colIdx = i * 128; tIdx < 128; ++tIdx, ++colIdx)
			{
				// generate the column indexed by colIdx. This is done with
				// AES in counter mode acting as a PRNG. We don't use the normal
				// PRNG interface because that would result in a data copy when
				// we move it into the T0,T1 matrices. Instead we do it directly.
				mGens[colIdx][0].ecbEncCounterMode(mGensBlkIdx[colIdx], superBlkSize, ((block*)t0.data() + superBlkSize * tIdx));
				mGens[colIdx][1].ecbEncCounterMode(mGensBlkIdx[colIdx], superBlkSize, ((block*)t1.data() + superBlkSize * tIdx));

				// increment the counter mode idx.
				mGensBlkIdx[colIdx] += superBlkSize;
			}

			// transpose our 128 columns of 1024 bits. We will have 1024 rows,
			// each 128 bits wide.
			transpose128x1024(t0);
			transpose128x1024(t1);

			// This is the index of where we will store the matrix long term.
			// doneIdx is the starting row. i is the offset into the blocks of 128 bits.
			// __restrict isn't crucial, it just tells the compiler that this pointer
			// is unique and it shouldn't worry about pointer aliasing.
			block* __restrict mT0Iter = mT0.data() + mT0.stride() * doneIdx + i;
			block* __restrict mT1Iter = mT1.data() + mT1.stride() * doneIdx + i;

			for (u64 rowIdx = doneIdx, j = 0; rowIdx < stopIdx; ++j)
			{
				// because we transposed 1024 rows, the indexing gets a bit weird. But this
				// is the location of the next row that we want. Keep in mind that we had long
				// **contiguous** columns.
				block* __restrict t0Iter = ((block*)t0.data()) + j;
				block* __restrict t1Iter = ((block*)t1.data()) + j;

				// do the copy!
				for (u64 k = 0; rowIdx < stopIdx && k < 128; ++rowIdx, ++k)
				{
					*mT0Iter = *(t0Iter);
					*mT1Iter = *(t1Iter);

					t0Iter += superBlkSize;
					t1Iter += superBlkSize;

					mT0Iter += mT0.stride();
					mT1Iter += mT0.stride();
				}
			}
		}

		doneIdx = stopIdx;
	}
	timer.setTimePoint("transpose T0 and T1");
	std::cout<<"finishing transposing T0 and T1"<<std::endl;


	int columnBytesNumber = (omega+7)/8;
	//send the correction to the sender
	std::vector<u8> sendBuffer(numOtExt*columnBytesNumber);

	u8* sendBuffer_ptr = sendBuffer.data();
	auto mT1_ptr = mT1.data();
	for(int i=0; i<numOtExt; i++){
		for(int j=0; j<columnBlockNumber; j++){
			mT1[i][j] ^= (mT0[i][j] ^ pax[i][j]);
		}
		std::memcpy(sendBuffer_ptr, (u8 *)(mT1_ptr), columnBytesNumber);
		sendBuffer_ptr += columnBytesNumber;
		mT1_ptr += mT1.stride();
	}
	cp::sync_wait(chl.send(sendBuffer));


	// std::cout<<"Receiver transpose the paxos table "<<std::endl;
	timer.setTimePoint("send the correction");
	double comSent_oprf = (chl.bytesSent())*1.0/(1<<20);
	double comRecv_oprf = (chl.bytesReceived())*1.0/(1<<20);
	lout << "Communication (MB): sent & recvd " << Color::Green << comSent_oprf << " " << comRecv_oprf << std::endl << Color::Default;
	std::cout<<"The OPRF communication (MB): "<<(comSent_oprf+comRecv_oprf)<<std::endl;

	//send the correction t0 ^ sMsgs[0] ^ sMsgs[1]
	oc::Matrix<block> vals(nx, columnBlockNumber);
	// paxos.template decode<block>(receiverSet, vals, mT0);
	paxos.decode<block>(receiverSet, vals, mT0, dt);

	// std::cout<<"Receiver transpose the paxos Q"<<std::endl;
	timer.setTimePoint("get the T0 values");

	// std::cout<<"The vals size: "<<vals.size()<<", "<<vals[0].size()<<std::endl;
	//computet the oprf values and put it into a map
	int outputByteNumber = valueBlocks*16;
	std::vector<std::vector<block>> receiverOPRFValues(nx, std::vector<block>(valueBlocks));

	RandomOracle H(outputByteNumber);
	for(int i=0; i<nx; i++){
		H.Reset();
		H.Update((u8 *)(vals.data()+i*vals.stride()), columnBytesNumber);
		H.Final((u8 *)receiverOPRFValues[i].data());
	}

	// std::cout<<"Receiver gets the oprf values"<<std::endl;

	timer.setTimePoint("get the oprf values");

	oc::Matrix<block> xorSums(paxos.size(), valueBlocks);
	cp::sync_wait(chl.recv(xorSums));

	oc::Matrix<block> devodedValues(nx, valueBlocks);
	// paxos.template decode<block>(receiverSet, devodedValues, xorSums);
	paxos.decode<block>(receiverSet, devodedValues, xorSums, dt);

	for(int i=0; i<nx; i++){
		for(int j=0; j<valueBlocks; j++){
			devodedValues[i][j] ^= receiverOPRFValues[i][j];
		}
	}

	timer.setTimePoint("get the opprf values");

	double comSent = (chl.bytesSent())*1.0/(1<<20);
	double comRecv = (chl.bytesReceived())*1.0/(1<<20);
	// std::cout<<"Receiving sender's oprf set (MB): "<<(comRecv - comRecv_oprf)<<std::endl;

	com = (chl.bytesSent() + chl.bytesReceived())*1.0/(1<<20);
	// lout << " n=" << Color::Green << omega << " " << milli << " ms  " << com << " MB" << std::endl << Color::Default;
	lout << "Communication (MB): sent & recvd " << Color::Green << comSent << " " << comRecv << std::endl << Color::Default;
	lout << "OPPRF communication (MB): " << (comSent + comRecv) << std::endl << Color::Default;


	lout << " **** receiver ****\n" << timer << std::endl;

	chl.close();
}

void runSenderOpprf(int nx, int ny, int valueBlocks, int omega, cp::AsioSocket &chl, bool verbose=0)
{
	oc::block setSeed1(toBlock(54321)), setSeed2(sysRandomSeed());
	std::vector<block> senderSet(ny);
	PRNG prng1_set(setSeed1), prng2_set(setSeed2);
	int intersectionSize = 100;
	prng1_set.get<block>(senderSet.data(), intersectionSize);
	prng2_set.get<block>(senderSet.data()+intersectionSize, ny-intersectionSize);

	std::vector<std::vector<block>> senderValues(ny, std::vector<block>(valueBlocks));
	for(int i=0; i<ny; i++){
		prng2_set.get<block>(senderValues[i].data(), valueBlocks);
	}

	Timer timer;

	auto s = timer.setTimePoint("start");

	int columnBlockNumber = (omega+127)/128;
	// std::vector<block> senderValues(n);

	//encode to to get okvs table;
	int w = 3;
	int ssp = 40;
	auto dt = PaxosParam::Binary;
	auto cols = 0;

	Baxos paxos;
	paxos.init(nx, 1<<15, w, ssp, dt, block(0, 0));

	// PaxosParam pp(nx, w, ssp, dt);
	std::cout<<"The expansion rate: "<<(paxos.size()*1.0/nx)<<std::endl;
	// oc::Matrix<block> valC(ny, columnBlockNumber), pax(paxos.size(), columnBlockNumber);

	oc::Matrix<block> valC(ny, columnBlockNumber);

	std::vector<block> valC_tmp(ny);
	std::vector<oc::AES> codewordsGenerator(columnBlockNumber);
	oc::block baseSeed(toBlock(11111));

	for(int i=0; i<codewordsGenerator.size(); i++){
		codewordsGenerator[i].setKey(baseSeed+toBlock(i));
		codewordsGenerator[i].ecbEncBlocks(senderSet, valC_tmp);
		for(int j=0; j<ny; j++){
			valC[j][i] = valC_tmp[j];
		}
	}

	PRNG prng(toBlock(12345));

	SilentOtExtReceiver  receiver;

	// Now compute the base OTs, we need to set them on the first pair of extenders.
	// In real code you would only have a sender or reciever, not both. But we do
	// here just showing the example.

	DefaultBaseOT base;
	std::array<std::array<block, 2>, 128> baseMsg;

	// perform the base To, call sync_wait to block until they have completed.
	cp::sync_wait(base.send(baseMsg, prng, chl));
	receiver.setBaseOts(baseMsg);

	u64 numCols = (omega+127)/128*128;

	BitVector choice(numCols);
	// in this case pick random messages.
	choice.randomize(prng);

	// construct a vector to stored the received messages.
	std::vector<block> rMsgs(numCols);
	// perform  totalOTs random OTs, the results will be written to msgs.
	cp::sync_wait(receiver.receive(choice, rMsgs, prng, chl));

	auto e = timer.setTimePoint("finish");
	auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(e - s).count();

	auto com = (chl.bytesSent() + chl.bytesReceived())*1.0/(1<<20);

	lout <<"ny="<<ny<< ", omega = " << Color::Green << numCols << " " << milli << " ms  " << com << " MB" << std::endl << Color::Default;

	int okvsTableSize = paxos.size();
	int rowBlockNumber = (okvsTableSize+127)/128;
	// std::vector<block> expand_msgs(rowBlockNumber);

	//the number of column should be divided by 128

	// round up
	int numOTExt = okvsTableSize;
	// numOTExt = ((numOTExt + 127) / 128) * 128;
	std::cout<<"The table size: "<<numOTExt<<std::endl;

	std::vector<AES> mGens(numCols);
	for(int i=0; i<numCols; i++){
		mGens[i].setKey(rMsgs[i]);
	}

	// We need two matrices, one for the senders matrix T^i_{b_i} and
	// one to hold the the correction values. This is sometimes called
	// the u = T0 + T1 + C matrix in the papers.
	oc::Matrix<block> mT(numOTExt, numCols/ 128);
	//char c;
	//chl.recv(&c, 1);

	oc::Matrix<block> mCorrectionVals(numOTExt, numCols / 128);

	// The receiver will send us correction values, this is the index of
	// the next one they will send.
	int mCorrectionIdx = 0;

	// we are going to process OTs in blocks of 128 * superblkSize messages.
	u64 numSuperBlocks = ((numOTExt + 127) / 128 + superBlkSize - 1) / superBlkSize;

	std::vector<u64> mGensBlkIdx(numCols, 0);

	// the index of the last OT that we have completed.
	u64 doneIdx = 0;

	// a temp that will be used to transpose the sender's matrix
	std::array<std::array<block, superBlkSize>, 128> t;

	for (u64 superBlkIdx = 0; superBlkIdx < numSuperBlocks; ++superBlkIdx)
	{
		// compute at what row does the user want use to stop.
		// the code will still compute the transpose for these
		// extra rows, but it is thrown away.
		u64 stopIdx
			= doneIdx
			+ std::min<u64>(u64(128) * superBlkSize, mT.bounds()[0] - doneIdx);

		// transpose 128 columns at at time. Each column will be 128 * superBlkSize = 1024 bits long.
		for (u64 i = 0; i < numCols / 128; ++i)
		{
			// generate the columns using AES-NI in counter mode.
			for (u64 tIdx = 0, colIdx = i * 128; tIdx < 128; ++tIdx, ++colIdx)
			{
				mGens[colIdx].ecbEncCounterMode(mGensBlkIdx[colIdx], superBlkSize, ((block*)t.data() + superBlkSize * tIdx));
				mGensBlkIdx[colIdx] += superBlkSize;
			}

			// transpose our 128 columns of 1024 bits. We will have 1024 rows,
			// each 128 bits wide.
			transpose128x1024(t);

			// This is the index of where we will store the matrix long term.
			// doneIdx is the starting row. i is the offset into the blocks of 128 bits.
			// __restrict isn't crucial, it just tells the compiler that this pointer
			// is unique and it shouldn't worry about pointer aliasing.
			block* __restrict mTIter = mT.data() + doneIdx * mT.stride() + i;

			for (u64 rowIdx = doneIdx, j = 0; rowIdx < stopIdx; ++j)
			{
				// because we transposed 1024 rows, the indexing gets a bit weird. But this
				// is the location of the next row that we want. Keep in mind that we had long
				// **contiguous** columns.
				block* __restrict tIter = (((block*)t.data()) + j);

				// do the copy!
				for (u64 k = 0; rowIdx < stopIdx && k < 128; ++rowIdx, ++k)
				{
					*mTIter = *tIter;

					tIter += superBlkSize;
					mTIter += mT.stride();
				}
			}

		}

		doneIdx = stopIdx;
	}

	timer.setTimePoint("transpose Q");
	std::cout<<"finishing transposing Q"<<std::endl;

	int columnByteNumber = (omega+7)/8;
	std::vector<u8> recvBuffer(columnByteNumber*numOTExt);
	cp::sync_wait(chl.recv(recvBuffer));

	std::vector<block> choiceBlocks(columnBlockNumber);
	std::memcpy((u8 *)choiceBlocks.data(), (u8 *)choice.data(), columnByteNumber);

	auto mCorrection_ptr = mCorrectionVals.data();
	auto recvBuffter_ptr = recvBuffer.data();
	for(int i=0; i<numOTExt; i++){
		// cp::sync_wait(chl.recv(recvBuffer));
		// std::memcpy((u8 *)(mCorrectionVals.data()+mCorrectionVals.stride()*i), recvBuffer.data()+columnByteNumber*i, columnByteNumber);
		std::memcpy((u8 *)(mCorrection_ptr), recvBuffter_ptr, columnByteNumber);
		mCorrection_ptr += mCorrectionVals.stride();
		recvBuffter_ptr += columnByteNumber;

		for(int j=0; j<columnBlockNumber; j++){
			mT[i][j] ^= (mCorrectionVals[i][j] & choiceBlocks[j]);
		}
	}
	oc::Matrix<block> vals1(ny, columnBlockNumber);
	// paxos.template decode<block>(senderSet, vals1, mT);
	paxos.decode<block>(senderSet, vals1, mT, dt);

	// std::cout<<"Sender gets the correction"<<std::endl;
	timer.setTimePoint("get correction");
	lout << "Communication (MB): sent & recvd " << Color::Green << (chl.bytesSent())*1.0/(1<<20) << " " << (chl.bytesReceived())*1.0/(1<<20) << std::endl << Color::Default;

	// int outputByteNumber = ((int)(ssp + log2(nx) + log2(ny))+7)/8;
	int outputByteNumber = valueBlocks*16;

	// std::cout<<"outputByteNumber: "<<outputByteNumber<<std::endl;
	RandomOracle H(outputByteNumber);
	std::vector<block> tmpResult(columnBlockNumber);
	oc::Matrix<block> xorSums(ny, valueBlocks), pax(paxos.size(), valueBlocks);

	for(int i=0; i<ny; i++){
		for(int j=0; j<columnBlockNumber; j++){
			tmpResult[j] = vals1[i][j] ^ (valC[i][j] & choiceBlocks[j]);
		}
		H.Reset();
		H.Update((u8 *)(tmpResult.data()), columnByteNumber);
		H.Final((u8 *)xorSums[i].data());
		for(int j=0; j<valueBlocks; j++){
			xorSums[i][j] ^= senderValues[i][j];
		}
	}
	timer.setTimePoint("get oprf values");

	// paxos.template encode<block>(xorSums, pax);
	paxos.solve<block>(senderSet, xorSums, pax, nullptr, dt);

	cp::sync_wait(chl.send(pax));

	timer.setTimePoint("sent opprf values");

	double comSent = (chl.bytesSent())*1.0/(1<<20);
	double comRecv = (chl.bytesReceived())*1.0/(1<<20);

	com = (chl.bytesSent() + chl.bytesReceived())*1.0/(1<<20);
	// lout << " n=" << Color::Green << omega << " " << milli << " ms  " << com << " MB" << std::endl << Color::Default;
	lout << "Communication (MB): sent & recvd " << Color::Green << comSent << " " << comRecv << std::endl << Color::Default;

	if (verbose)
	{
		lout << " **** sender ****\n" << timer << std::endl;
	};
	chl.close();
}

int main(int argc, char** argv)
{
    oc::CLP cmd(argc, argv);

	Role role = Role(cmd.get<int>("role"));

	// auto nr = cmd.getOr("receiverSize", 100);
	auto verbose = cmd.isSet("v");
	
	//the set size of the receiver and the sender are respectively: nx, ny
	int nx = 65536, ny = 1<<16;
	if(cmd.isSet("nnx")){
		nx = 1<<cmd.get<int>("nnx");
	}else if(cmd.isSet("nx")){
		nx = cmd.get<int>("nx");
	}
	if(cmd.isSet("nny")){
		ny = 1<<cmd.get<int>("nny");
	}else if(cmd.isSet("ny")){
		ny = cmd.get<int>("ny");
	}

	// bool opprf = cmd.isSet("opprf");
	int valueBlocks = 0;
	if(cmd.isSet("theta")){
		valueBlocks = cmd.get<int>("theta");
	}
	
	int omega = 444;
	if(nx == (1<<12)){
		omega = 420;   
	}else if(nx == (1<<16)){
		omega = 428;
	}else if(nx == (1<<20)){
		omega = 436;
	}else if(nx == (1<<24)){
		omega = 444;
	}
	
	// The statistical security parameter.
	// auto ssp = cmd.getOr("ssp", 40);

	std::string ip = "127.0.0.1:12121";

	if(role == Role::Receiver){
		// get up the networking, 1 for the receiver
		auto chl = cp::asioConnect(ip, 1);
		// runReceiver(nx, ny, omega, chl, verbose);
		if(valueBlocks > 0){
			std::cout<<"In receiver OPPRF: "<<std::endl;
			runReceiverOpprf(nx, ny, valueBlocks, omega, chl, verbose);
		}else{
			runReceiver(nx, ny, omega, chl, verbose);
		}
	}else{
		auto chl = cp::asioConnect(ip, 0);
		// runSender(nx, ny, omega, chl, verbose);
		if(valueBlocks > 0){
			std::cout<<"In sender OPPRF: "<<std::endl;
			runSenderOpprf(nx, ny, valueBlocks, omega, chl, verbose);
		}else{
			runSender(nx, ny, omega, chl, verbose);
		}
	}
	// run(role, totalOTs, ip, cmd);

	return 0;
}