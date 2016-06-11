#include "movePicker.hpp"
#include "generateMoves.hpp"
#include "thread.hpp"

namespace {
const Score LVATable[PieceTypeNum] = {
	Score(0), Score(1), Score(2), Score(3), Score(4), Score(7), Score(8), Score(6), Score(10000),
	Score(5), Score(5), Score(5), Score(5), Score(9), Score(10)
};
inline Score LVA(const PieceType pt) { return LVATable[pt]; }

  void insertion_sort(MoveStack* begin, MoveStack* end)
  {
    MoveStack tmp, *p, *q;

    for (p = begin + 1; p < end; ++p)
    {
        tmp = *p;
        for (q = p; q != begin && *(q-1) < tmp; --q)
            *q = *(q-1);
        *q = tmp;
    }
  }

  inline Move pick_best(MoveStack* begin, MoveStack* end)
  {
      std::swap(*begin, *std::max_element(begin, end));
      return *begin;
  }

} // namespace

MovePicker::MovePicker(const Position& p, const Move ttm, const Depth d, Search::Stack* s)
	: pos(p), ss(s), depth(d)
{
	assert(Depth0 < d);

	endBadCaptures = moves + MaxLegalMoves - 1;
    Square prevSq = (ss-1)->currentMove.to();
    countermove = pos.thisThread()->counterMoves[pos.piece(prevSq)][prevSq];

    stage = pos.inCheck() ? EvasionSearch : MainSearch;

	ttMove = (!ttm.isNone() && pos.moveIsPseudoLegal(ttm) ? ttm : Move::moveNone());
    endMoves += (!ttMove.isNone());
}

// 静止探索で呼ばれる。
MovePicker::MovePicker(const Position& p, Move ttm, const Depth d, const Square sq)
	: pos(p)
{
	assert(d <= Depth0);

	if (pos.inCheck())
        stage = QEvasionSearch;

	// todo: ここで Stockfish は qcheck がある。

	else if (DepthQRecaptures < d)
        stage = QSearch;

	else {
        stage = QRecapture;
		recaptureSquare = sq;
		ttm = Move::moveNone();
	}

	ttMove = (!ttm.isNone() && pos.moveIsPseudoLegal(ttm) ? ttm : Move::moveNone());
    endMoves += !ttMove.isNone();
}

MovePicker::MovePicker(const Position& p, const Move ttm, Score th)
	: pos(p), threshold(th)
{
	assert(!pos.inCheck());

    stage = ProbCut;

    ttMove = !ttm.isNone()
      && pos.moveIsPseudoLegal(ttm)
      && !ttMove.isCapture()
      && pos.see(ttm) > threshold ? ttm : MOVE_NONE;

    endMoves += !ttMove.isNone();
}

void MovePicker::scoreCaptures() {
	for (auto& m : *this) {
		const Move move = m;
		m.score = Position::pieceScore(pos.piece(move.to())) - LVA(move.pieceTypeFrom());
	}
}

template <bool IsDrop> void MovePicker::scoreNonCapturesMinusPro() {
  const HistoryStats& history = pos.thisThread()->history;

  const CounterMoveStats* cm = (ss-1)->counterMoves;
  const CounterMoveStats* fm = (ss-2)->counterMoves;
  const CounterMoveStats* f2 = (ss-4)->counterMoves;

	for (auto& m : *this) {
		const Move move = m;

        m.score = history[pos.moved_piece(move)][move.to()]
          + (cm ? (*cm)[pos.moved_piece(move)][move.to()] : ScoreZero)
          + (fm ? (*fm)[pos.moved_piece(move)][move.to()] : ScoreZero)
          + (f2 ? (*f2)[pos.moved_piece(move)][move.to()] : ScoreZero);
	}
}

void MovePicker::scoreEvasions() {
  const HistoryStats& history = pos.thisThread()->history;

	for (auto& m : *this) {
		const Move move = m;
		const Score seeScore = pos.seeSign(move);
		if (seeScore < 0)
			m.score = seeScore - HistoryStats::Max;
		else if (move.isCaptureOrPromotion()) {
			m.score = pos.capturePieceScore(pos.piece(move.to())) + HistoryStats::Max;
			if (move.isPromotion()) {
				const PieceType pt = pieceToPieceType(pos.piece(move.from()));
				m.score += pos.promotePieceScore(pt);
			}
		}
		else
			//it->score = history.value(move.isDrop(), colorAndPieceTypeToPiece(pos.turn(), move.pieceTypeFromOrDropped()), move.to());
            m.score = history.value(move.isDrop(), pos.moved_piece(move), move.to());
	}
}

void MovePicker::generate_next_stage() {

    assert(stage != PH_Stop);

    cur = moves;

	switch (++stage) {

	case PH_TacticalMoves0: case PH_TacticalMoves1:
        endMoves = generateMoves<CapturePlusPro>(cur, pos);
		scoreCaptures();
        break;

	case PH_Killers:
        killerMoves[0] = ss->killers[0];
        killerMoves[1] = ss->killers[1];
        killerMoves[2] = countermove;
        cur = killerMoves;
        endMoves = cur + 2 + (countermove != killerMoves[0] && countermove != killerMoves[1]);
        break;

	case PH_NonTacticalMoves:
        endMoves = generateMoves<NonCaptureMinusPro>(cur, pos);
		scoreNonCapturesMinusPro<false>();
		cur = endMoves;
		endMoves = generateMoves<Drop>(cur, pos);
		scoreNonCapturesMinusPro<true>();
		cur = moves;
        if (depth < static_cast<Depth>(3 * OnePly))
        {
            MoveStack* goodQuiet = std::partition(cur, endMoves, [](const MoveStack& m)
                                                 { return m.score > ScoreZero; });
            insertion_sort(cur, goodQuiet);
        }
        else
            insertion_sort(cur, endMoves);
        break;

	case PH_BadCaptures:
		cur = moves + MaxLegalMoves - 1;
        endMoves = endBadCaptures;
        break;

	case PH_Evasions: case PH_QEvasions:
        endMoves = generateMoves<Evasion>(cur, pos);
		if (endMoves - moves > 1)
			scoreEvasions();
        break;

	case PH_QCaptures0:
        endMoves = generateMoves<CapturePlusPro>(moves, pos);
		scoreCaptures();
        break;

	case PH_QCaptures1:
        endMoves = generateMoves<Recapture>(moves, pos, recaptureSquare);
		scoreCaptures();
        break;

    case EvasionSearch: case QSearch: case QEvasionSearch: case QRecapture: case ProbCut: case PH_Stop:
		// これが無いと、MainSearch の後に EvasionSearch が始まったりしてしまう。
        stage = PH_Stop;
        break;

	default: UNREACHABLE;
	}
}

Move MovePicker::nextMove() {
	Move move;
	
    while (true)
    {
		// end() に達したら次の phase に移る。
		while (cur == endMoves && stage != PH_Stop)
            generate_next_stage();

		switch (stage) {

		case MainSearch: case EvasionSearch: case QSearch: case QEvasionSearch: case ProbCut:
			++cur;
			return ttMove;

		case PH_TacticalMoves0:
			move = pick_best(cur++, endMoves);
			if (move != ttMove) {
				if (ScoreZero <= pos.see(move))
					return move;

				// 後ろから SEE の点数が高い順に並ぶようにする。
				*endBadCaptures-- = move;
			}
			break;

		case PH_Killers:
			move = *cur++;
			if (!move.isNone()
				&& move != ttMove
				&& pos.moveIsPseudoLegal(move, true)
				&& pos.piece(move.to()) == Empty)
				return move;
			break;

		case PH_NonTacticalMoves:
			move = *cur++;
			if (move != ttMove
				&& move != killerMoves[0]
                && move != killerMoves[1]
                && move != killerMoves[2])
				return move;
			break;

		case PH_BadCaptures:
			return *cur--;

		case PH_Evasions: case PH_QEvasions: case PH_QCaptures0:
			move = pick_best(cur++, endMoves);
			if (move != ttMove)
				return move;
			break;

		case PH_TacticalMoves1:
			move = pick_best(cur++, endMoves);
			// todo: see が確実に駒打ちじゃないから、内部で駒打ちか判定してるのは少し無駄。
			if (move != ttMove && threshold < pos.see(move))
				return move;
			break;

		case PH_QCaptures1:
			move = pick_best(cur++, endMoves);
			assert(move.to() == recaptureSquare);
			return move;

		case PH_Stop:
			return Move::moveNone();

		default:
			UNREACHABLE;
		}
	}
}
